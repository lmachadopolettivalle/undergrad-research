import argparse

import tangos as db
import tangos.core
from tangos import parallel_tasks
from tangos import core
import sqlalchemy, sqlalchemy.orm
from tangos.log import logger
import numpy as np


class GenericLinker(object):
    def __init__(self, session=None):
        self.session = session or core.get_default_session()

    def parse_command_line(self, argv=None):
        parser = self._create_argparser()
        self.args = parser.parse_args(argv)
        db.process_options(self.args)

    def _create_argparser(self):
        parser = argparse.ArgumentParser()
        db.supplement_argparser(parser)
        parser.add_argument("--verbose", action="store_true",
                            help="Print extra information")
        parser.add_argument("--force", action="store_true",
                            help="Generate links even if they already exist for those timesteps")
        parser.add_argument("--hmax", action="store", type=int, default=None,
                            help="Specify the maximum number of halos per snapshot")
        parser.add_argument('--backwards', action='store_true',
                            help='Process in reverse order (low-z first)')
        parser.add_argument('--dmonly', action='store_true',
                            help='only match halos based on DM particles. Much more memory efficient, but currently only works for Rockstar halos')
        return parser

    def run_calculation_loop(self):

        pair_list = self._generate_timestep_pairs()

        if len(pair_list)==0:
            logger.error("No timesteps found to link")
            return

        pair_list = parallel_tasks.distributed(pair_list)

        for s_x, s in pair_list:
            logger.info("Linking %r and %r",s_x,s)
            if self.args.force or self.need_crosslink_ts(s_x, s):
                self.crosslink_ts(s_x, s, 0, self.args.hmax, self.args.dmonly)

    def _generate_timestep_pairs(self):
        raise NotImplementedError, "No implementation found for generating the timestep pairs"

    def get_halo_entry(self, ts, halo_number):
        h = ts.halos.filter_by(finder_id=halo_number).first()
        return h

    def need_crosslink_ts(self, ts1, ts2):
        num_sources = ts1.halos.count()
        num_targets = ts2.halos.count()
        if num_targets == 0:
            logger.warn("Will not link: no halos in target timestep %r", ts2)
            return False
        if num_sources == 0:
            logger.warn("Will not link: no halos in source timestep %r", ts1)
            return False

        halo_source = sqlalchemy.orm.aliased(core.halo.Halo, name="halo_source")
        halo_target = sqlalchemy.orm.aliased(core.halo.Halo, name="halo_target")
        same_d_id = core.dictionary.get_or_create_dictionary_item(self.session, "ptcls_in_common").id
        exists = self.session.query(core.halo_data.HaloLink).join(halo_source, core.halo_data.HaloLink.halo_from). \
                    join(halo_target, core.halo_data.HaloLink.halo_to). \
                    filter(halo_source.timestep_id == ts1.id, halo_target.timestep_id == ts2.id,
                        core.halo_data.HaloLink.relation_id == same_d_id).count() > 0
        self.session.commit()

        if exists:
            logger.warn("Will not link: links already exist between %r and %r", ts1, ts2)
            return False
        return True

    def create_db_objects_from_catalog(self, cat, halos1, halos2, fid1, fid2, same_d_id):
        items = []
        missing_db_object = 0
        for i, possibilities in enumerate(cat):
            o1 = np.where(fid1==i)[0]
            if len(o1)>0:
                h1 = halos1[o1[0]]
            else:
                h1 = None
            for cat_i, weight in possibilities:
                o2 = np.where(fid2==cat_i)[0]
                if len(o2)>0:
                    h2 = halos2[o2[0]]
                else:
                    h2 = None
                if h1 is not None and h2 is not None:
                    items.append(core.halo_data.HaloLink(h1, h2, same_d_id, weight))
                else:
                    missing_db_object += 1

        if missing_db_object > 0:
            logger.warn("%d link(s) could not be identified because the halo objects do not exist in the DB",
                        missing_db_object)
        return items

    def crosslink_ts(self, ts1, ts2, halo_min=0, halo_max=None, dmonly=False, threshold=0.005):
        """Link the halos of two timesteps together

        :type ts1 tangos.core.TimeStep
        :type ts2 tangos.core.TimeStep"""
        logger.info("Gathering halo information for %r and %r", ts1, ts2)
        halos1 = ts1.halos.all()
        halos2 = ts2.halos.all()
        fid1 = np.array([h.finder_id for h in halos1])
        fid2 = np.array([h.finder_id for h in halos2])

        with parallel_tasks.RLock("create_db_objects_from_catalog"):
            same_d_id = core.dictionary.get_or_create_dictionary_item(self.session, "ptcls_in_common")
            self.session.commit()

        output_handler_1 = ts1.simulation.get_output_set_handler()
        output_handler_2 = ts2.simulation.get_output_set_handler()
        if not isinstance(output_handler_1, type(output_handler_2)):
            logger.error("Timesteps %r and %r cannot be crosslinked; they are using incompatible file readers",
                         ts1, ts2)
            return

        snap1 = ts1.load()
        snap2 = ts2.load()

        try:
            cat = output_handler_1.match_halos(snap1, snap2, halo_min, halo_max, dmonly, threshold)
            back_cat = output_handler_2.match_halos(snap2, snap1, halo_min, halo_max, dmonly, threshold)
        except:
            logger.exception("Exception during attempt to crosslink timesteps %r and %r", ts1, ts2)
            return

        with self.session.no_autoflush:
            logger.info("Gathering links for %r and %r", ts1, ts2)
            items = self.create_db_objects_from_catalog(cat, halos1, halos2, fid1, fid2, same_d_id)
            logger.info("Identified %d links between %r and %r", len(items), ts1, ts2)
            items_back = self.create_db_objects_from_catalog(back_cat, halos2, halos1, fid2, fid1, same_d_id)
            logger.info("Identified %d links between %r and %r", len(items_back), ts2, ts1)

        with parallel_tasks.RLock("create_db_objects_from_catalog"):
            logger.info("Preparing to commit links for %r and %r", ts1, ts2)
            self.session.add_all(items)
            self.session.add_all(items_back)
            self.session.commit()
        logger.info("Finished committing total of %d links for %r and %r", len(items)+len(items_back), ts1, ts2)

class TimeLinker(GenericLinker):
    def _generate_timestep_pairs(self):
        logger.info("generating pairs of timesteps")
        base_sim = db.sim_query_from_name_list(self.args.sims)
        pairs = []
        for x in base_sim:
            ts = self.session.query(tangos.core.timestep.TimeStep).filter_by(
                simulation_id=x.id, available=True).order_by(tangos.core.timestep.TimeStep.redshift.desc()).all()
            for a, b in zip(ts[:-1], ts[1:]):
                pairs.append((a, b))
        if self.args.backwards:
            pairs = pairs[::-1]
        return pairs

    def _create_argparser(self):
        parser = super(TimeLinker, self)._create_argparser()
        parser.add_argument('--sims', '--for', action='store', nargs='*',
                            metavar='simulation_name',
                            help='Specify a simulation (or multiple simulations) to run on')
        return parser

class CrossLinker(GenericLinker):
    def _create_argparser(self):
        parser = super(CrossLinker, self)._create_argparser()
        parser.add_argument('sims', action='store', nargs=2,
                            metavar=('name1', 'name2'),
                            help='The two simulations (or, optionally, two individual timesteps) to crosslink')
        return parser

    def _generate_timestep_pairs(self):
        obj1 = db.get_item(self.args.sims[0])
        obj2 = db.get_item(self.args.sims[1])
        if isinstance(obj1, core.TimeStep) and isinstance(obj2, core.TimeStep):
            return [[obj1, obj2]]
        elif isinstance(obj1, core.Simulation) and isinstance(obj2, core.Simulation):
            return self._generate_timestep_pairs_from_sims(obj1, obj2)
        else:
            raise ValueError, "No way to link %r to %r"%(obj1, obj2)

    def _generate_timestep_pairs_from_sims(self, sim1, sim2):
        assert sim1 != sim2, "Can't link simulation to itself"

        logger.info("Match timesteps of %r to %r",sim1,sim2)

        ts1s = sim1.timesteps
        ts2s = sim2.timesteps

        pairs = []
        for ts1 in ts1s:
            ts2 = self._get_best_timestep_matching(ts2s, ts1)
            pairing_is_mutual = (self._get_best_timestep_matching(ts1s,ts2)==ts1)
            if pairing_is_mutual:
                logger.info("Pairing timesteps: %r and %r", ts1, ts2)
                pairs+=[(ts1,ts2)]
            else:
                logger.warn("No pairing found for timestep %r",ts1)

        return pairs

    def _get_best_timestep_matching(self, list_of_candidates, timestep_to_match):
        return min(list_of_candidates, key=lambda ts2: abs(ts2.time_gyr - timestep_to_match.time_gyr))

import weakref

from sqlalchemy import Column, Integer, ForeignKey, orm
from sqlalchemy.orm import relationship, backref

from . import extraction_patterns
from . import Base
from . import creator
from .dictionary import get_dict_id, get_or_create_dictionary_item
from .timestep import TimeStep
from .tracking import TrackerHaloCatalogue


class Halo(Base):
    __tablename__ = 'halos'


    id = Column(Integer, primary_key=True)
    halo_number = Column(Integer)
    finder_id = Column(Integer)
    timestep_id = Column(Integer, ForeignKey('timesteps.id'))
    timestep = relationship(TimeStep, backref=backref(
        'halos', order_by=halo_number, cascade='all', lazy='dynamic'), cascade='save-update, merge')
    NDM = Column(Integer)
    NStar = Column(Integer)
    NGas = Column(Integer)
    creator = relationship(creator.Creator, backref=backref(
        'halos', cascade='all', lazy='dynamic'), cascade='save-update')
    creator_id = Column(Integer, ForeignKey('creators.id'))
    halo_type = Column(Integer, nullable=False)

    __mapper_args__ = {
        'polymorphic_identity':0,
        'polymorphic_on':halo_type
    }

    def __init__(self, timestep, halo_number, finder_id, NDM, NStar, NGas, halo_type=0):
        self.timestep = timestep
        self.halo_number = halo_number
        self.finder_id = finder_id
        self.NDM = NDM
        self.NStar = NStar
        self.NGas = NGas
        self.halo_type = halo_type
        self.init_on_load()
        self.creator_id = creator.get_creator_id()

    @orm.reconstructor
    def init_on_load(self):
        self._dict_is_complete = False
        self._d = {}

    def __repr__(self):

        return "<%s %r | NDM=%d Nstar=%d Ngas=%d>"%(self.__class__.__name__, self.path, self.NDM, self. NStar, self.NGas)


    @property
    def path(self):
        if self.halo_type==0:
            name = str(self.halo_number)
        else:
            name = str(self.halo_type)+"."+str(self.halo_number)
        return self.timestep.path+"/"+name


    def load(self, partial=False):
        """Use pynbody to load the data for this halo, if it is present on this computer's filesystem.

        By default, the entire simulation is loaded and a subview with this halo is returned. If partial=True,
        pynbody's partial loading system is used to only load the data for the one halo, saving memory."""

        handler = self.timestep.simulation.get_output_set_handler()
        if self.halo_type==0:
            return handler.load_halo(self.timestep.extension, self.finder_id, partial=partial)
        else:
            return handler.load_tracked_region(self.timestep.extension, self.tracker, partial=partial)


    def calculate(self, name):
        """Use the live-calculation system to calculate a user-specified function of the stored data.

        See live_calculation.md for an introduction to this powerful functionality."""
        from .. import live_calculation
        calculation = live_calculation.parser.parse_property_name_if_required(name)
        value = calculation.values_sanitized([self])[0]
        if len(value)==1:
            return value[0]
        else:
            return value

    def calculate_abcissa_values(self, name):
        from .. import live_calculation
        calculation = live_calculation.parser.parse_property_name_if_required(name)
        (value, ), description = calculation.values_and_description([self])
        return description.plot_x_values(value[0])

    def __getitem__(self, key):
        """Highest-level method for retrieving data or link"""
        return self.get_data(key)

    def get(self, key, default=None):
        """Highest-level method for retrieving data or link, with default return value if no data is found
        under the specified key"""
        try:
            return self.get_data(key)
        except KeyError:
            return default

    def get_data(self, key, raw=False, always_return_array=False):
        """High-level data access to the property or linked halo named by key.

        :param key: string with the name of the property or link
        :param raw: if True, get the raw data rather than attempt to reassemble the data into a science-ready form
        :param always_return_array: if True, always return the data as a list of values, even if there is only one
           stored property with the specified name
        """

        if raw:
            getters=[extraction_patterns.halo_property_raw_value_getter]
        else:
            getters=[extraction_patterns.halo_property_value_getter]
        getters+=[extraction_patterns.halo_link_target_getter]

        return_data = self.get_objects(key, getters)

        if (not always_return_array) and len(return_data) == 1:
            return_data = return_data[0]

        return return_data

    def get_objects(self, key, getters = [extraction_patterns.halo_property_getter,
                                          extraction_patterns.halo_link_getter]):
        """Get objects belonging to this halo named by the specified key.

        Compared to get_data, this allows access to the underlying HaloProperty or HaloLink objects, or to perform
        custom processing of the data by specifying particular extraction patterns to getters. For more information,
        see halo_data_extraction_patterns."""
        from . import Session
        session = Session.object_session(self)
        key_id = get_dict_id(key, session=session)

        ret_values = []
        for g in getters:
            ret_values += g.get(self, key_id, session)

        if len(ret_values) == 0:
            raise KeyError, "No such property %r" % key
        return ret_values


    def __setitem__(self, key, obj):
        if isinstance(obj, Halo):
            self._setitem_one_halo(key, obj)
        elif hasattr(obj, '__len__') and all([isinstance(x,Halo) for x in obj]):
            self._setitem_multiple_halos(key, obj)
        else:
            self._setitem_property(key, obj)

    def _setitem_property(self, key, obj):
        from . import Session
        from .halo_data import HaloProperty

        session = Session.object_session(self)
        key = get_or_create_dictionary_item(session, key)
        X = self.properties.filter_by(name_id=key.id).first()
        if X is not None:
            X.data = obj
        else:
            X = session.merge(HaloProperty(self, key, obj))
        X.creator_id = creator.get_creator_id()

    def _setitem_one_halo(self, key, obj):
        from . import Session, HaloLink
        session = Session.object_session(self)
        key = get_or_create_dictionary_item(session, key)
        X = self.links.filter_by(halo_from_id=self.id, relation_id=key.id).first()
        if X is None:
            X = session.merge(HaloLink(self, obj, key))
            X.creator_id = creator.get_creator_id()
        else:
            X.halo_to = obj

    def _setitem_multiple_halos(self, key, obj):
        from . import Session
        from .halo_data import HaloLink

        session = Session.object_session(self)
        key = get_or_create_dictionary_item(session, key)
        self.links.filter_by(halo_from_id=self.id, relation_id=key.id).delete()
        links = [HaloLink(self, halo_to, key) for halo_to in obj]
        session.add_all(links)


    def keys(self, getters = [extraction_patterns.halo_property_getter,
                              extraction_patterns.halo_link_getter]):
        from . import Session
        names = []
        session = Session.object_session(self)
        for g in getters:
            names+=g.keys(self, session)


        return names

    def __contains__(self, item):
        return item in self.keys()

    @property
    def tracker(self):
        if self.halo_type != 1:
            return None
        else:
            return self.timestep.simulation.trackers.filter_by(halo_number=self.halo_number).first()

    @property
    def earliest(self):
        if self.previous is not None:
            return self.previous.earliest
        else:
            return self

    @property
    def latest(self):
        if self.next is not None:
            return self.next.latest
        else:
            return self

    def plot(self, name, *args, **kwargs):
        from . import Session
        name_id = get_dict_id(name, Session.object_session(self))
        data = self.properties.filter_by(name_id=name_id).first()
        return data.plot(*args, **kwargs)


    def property_cascade(self, *plist, **kwargs):
        """Run the specified calculations on this halo and its descendants

        Each argument is a string (or an instance of live_calculation.Calculation), following the syntax
        described in live_calculation.md.

        *kwargs*:

        :param nmax: The maximum number of descendants to consider (default 1000)
        :param strategy: The class to use to find the descendants (default relation_finding.MultiHopMajorDescendantsStrategy)
        """
        from .. import live_calculation
        from .. import relation_finding
        from .. import temporary_halolist as thl
        from . import Session
        from .. import query as db_query

        nmax = kwargs.get('nmax',1000)
        strategy = kwargs.get('strategy', relation_finding.MultiHopMajorDescendantsStrategy)
        strategy_kwargs = kwargs.get('strategy_kwargs', {})

        if isinstance(plist[0], live_calculation.Calculation):
            property_description = plist[0]
        else:
            property_description = live_calculation.parser.parse_property_names(*plist)

        # must be performed in its own session as we intentionally load in a lot of
        # objects with incomplete lazy-loaded properties
        session = Session()
        try:
            with strategy(db_query.get_halo(self.id, session), nhops_max=nmax,
                          include_startpoint=True, **strategy_kwargs).temp_table() as tt:
                raw_query = thl.halo_query(tt)
                query = property_description.supplement_halo_query(raw_query)
                results = query.all()
                return property_description.values_sanitized(results)
        finally:
            session.close()

    def reverse_property_cascade(self, *plist, **kwargs):
        """Run the specified calculations on the progenitors of this halo

        For more information see property_cascade.
        """
        from .. import relation_finding
        kwargs['strategy'] = relation_finding.MultiHopMajorProgenitorsStrategy
        return self.property_cascade(*plist, **kwargs)


    @property
    def next(self):
        if not hasattr(self, '_next'):
            from .. import relation_finding
            strategy = relation_finding.HopMajorDescendantStrategy(self)
            self._next=strategy.first()

        return self._next

    @property
    def previous(self):
        if not hasattr(self, '_previous'):
            from .. import relation_finding
            strategy = relation_finding.HopMajorProgenitorStrategy(self)
            self._previous=strategy.first()

        return self._previous

    def short(self):
        return "<Halo " + str(self.halo_number) + " of ...>"


class BH(Halo):
    __mapper_args__ = {
        'polymorphic_identity':1
    }

    def __init__(self, timestep, halo_number):
        super(BH, self).__init__(timestep, halo_number, halo_number, 0,0,0,1)


_loaded_halocats = {}

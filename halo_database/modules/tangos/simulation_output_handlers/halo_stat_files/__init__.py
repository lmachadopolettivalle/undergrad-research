import os

import numpy as np

import tangos.core.dictionary
import tangos.core.halo
import tangos.core.halo_data
from tangos import core
from . import translations


class HaloStatFile(object):
    """Manages and reads a halo stat file of unspecified format."""
    _id_offset = 0
    _column_translations = {}

    def __new__(cls, timestep):
        for subclass in cls.__subclasses__():
            if subclass.can_load(timestep):
                return object.__new__(subclass, timestep)
        else:
            raise IOError, "No stat file found for this timestep"

    @classmethod
    def can_load(cls, timestep):
        return os.path.exists(cls.filename(timestep))

    @classmethod
    def filename(cls, timestep):
        raise ValueError, "Unknown path to stat file"

    def __init__(self, timestep):
        self._timestep = timestep
        self.filename = self.filename(timestep)

    def iter_rows_raw(self, *args):
        """
        Yield the halo ID and requested columns from each line of the stat file, without any emulation.

        :param args: strings for the column names
        :return: id, arg1, arg2, arg3 where ID is the halo ID and argN is the value of the Nth named column
        """
        with open(self.filename) as f:
            header = self._read_column_names(f)
            ids = [0] + [header.index(a) for a in args]
            for l in f:
                yield self._get_values_for_columns(ids, l)

    def iter_rows(self, *args):
        """
        Yield the halo ID and requested columns from each line of the stat file, emulating the existence of some columns.

        For example, AHF stat files do not contain n_dm; however, its value can be automatically inferred by this
        function. Meanwhile IDL .amiga.stat files rename n_gas as N_gas.

        :param args: strings for the column names
        :return: id, arg1, arg2, arg3 where ID is the halo ID and argN is the value of the Nth named column
        """

        raw_args = []
        for arg in args:
            if arg in self._column_translations:
                raw_args+=self._column_translations[arg].inputs()
            else:
                raw_args.append(arg)

        for raw_values in self.iter_rows_raw(*raw_args):
            values = [raw_values[0]]
            for arg in args:
                if arg in self._column_translations:
                    values.append(self._column_translations[arg](raw_args, raw_values[1:]))
                else:
                    values.append(raw_values[1:][raw_args.index(arg)])
            yield values

    def read(self, *args):
        """Read the halo ID and requested columns from the entire file, returning each column as a separate array"""
        return_values = [[] for _ in xrange(len(args)+1)]
        for row in self.iter_rows(*args):
            for return_array, value in zip(return_values, row):
                return_array.append(value)

        return [np.array(x) for x in return_values]

    def _get_values_for_columns(self, columns, line):
        results = []
        l_split = line.split()
        for id_this in columns:
            this_str = l_split[id_this]
            if "." in this_str:
                guess_type = float
            else:
                guess_type = int

            try:
                this_cast = guess_type(this_str)
            except ValueError:
                this_cast = this_str

            results.append(this_cast)
        results[0] += self._id_offset
        return results

    def _read_column_names(self, f):
        return [x.split("(")[0] for x in f.readline().split()]

    def make_halo_objects(self, min_NDM=1000):
        halos = []
        import numpy as np
        n_tot = []
        for num, NDM, Nstar, Ngas in self.iter_rows("n_dm", "n_star", "n_gas"):
            n_tot.append(NDM+Nstar+Ngas)
        ids = np.zeros(len(n_tot))
        ids[np.argsort(np.array(n_tot))[::-1]] = np.arange(len(n_tot)) + 1
        cnt = 1
        for num, NDM, Nstar, Ngas in self.iter_rows("n_dm", "n_star", "n_gas"):
            if NDM > min_NDM:
                h = tangos.core.halo.Halo(self._timestep, ids[cnt-1],cnt, NDM, Nstar, Ngas)
                halos.append(h)
            cnt += 1
        return halos

    def add_halos(self, min_NDM=1000):
        halos = self.make_halo_objects(min_NDM)
        session = core.Session.object_session(self._timestep)
        session.add_all(halos)

    def add_halo_properties(self, *property_names):
        halos = self._timestep.halos.all()
        halos_map = {}
        for h in halos:
            halos_map[h.halo_number] = h

        session = core.Session.object_session(self._timestep)

        property_db_names = [tangos.core.dictionary.get_or_create_dictionary_item(session, name) for name in property_names]
        property_objects = []
        for values in self.iter_rows(*property_names):
            halo = halos_map.get(values[0], None)
            if halo is not None:
                for name_object, value in zip(property_db_names, values[1:]):
                    property_objects.append(tangos.core.halo_data.HaloProperty(halo, name_object, value))

        session.add_all(property_objects)


class AHFStatFile(HaloStatFile):
    _id_offset = 1

    _column_translations = {'n_dm': translations.Function(lambda ngas, nstar, npart: npart - ngas - nstar,
                                                          'n_gas', 'n_star', 'npart')}

    @classmethod
    def filename(cls, timestep):
        import glob
        file = glob.glob(timestep.filename+'.z*.???.AHF_halos')[0]
        return file


class AmigaIDLStatFile(HaloStatFile):
    _id_offset = 0

    _column_translations = {'n_dm': translations.Rename('N_dark'),
                            'n_gas': translations.Rename('N_gas'),
                            'n_star': translations.Rename("N_star"),
                            'npart': translations.Function(lambda ngas, nstar, ndark: ngas + nstar + ndark,
                                                           "N_dark", "N_gas", "N_star")}

    @classmethod
    def filename(cls, timestep):
        return timestep.filename + '.amiga.stat'


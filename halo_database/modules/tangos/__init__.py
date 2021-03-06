"""tangos - the amazing numerical galaxy organization system

This package provides ways to create and interpret databases of numerical galaxy formation simulations for scientific
analysis.

For information on getting started, see README.md.

"""

import sqlalchemy
import sqlalchemy.orm.session
from sqlalchemy import Index, Column, Integer, String, Float, ForeignKey, DateTime, Boolean, LargeBinary, create_engine, orm
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref, sessionmaker, clear_mappers, deferred
from sqlalchemy import and_, or_
from sqlalchemy.orm.session import Session




from . import core, log
from .core import *
from .query import *

core.init_db()
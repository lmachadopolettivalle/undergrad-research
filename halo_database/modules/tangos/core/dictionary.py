import sqlalchemy, sqlalchemy.exc
from sqlalchemy import Column, Integer, String

from . import Base, get_default_session

_dict_id = {}  # maps dictionary text -> database ID
_dict_obj = {} # maps session, dictionary text -> database object


class DictionaryItem(Base):
    __tablename__ = 'dictionary'

    id = Column(Integer, primary_key=True)
    text = Column(String, unique=True)

    def __repr__(self):
        return "<DictionaryItem " + self.text + ">"

    def __init__(self, text):
        self.text = text

    def providing_class(self):
        import properties
        return properties.providing_class(self.text)

def get_dict_id(text, default=None, session=None):
    """Get a DictionaryItem id for text (possibly cached). Raises KeyError if
    no dictionary object exists for the specified text"""

    from . import Session

    if session is None:
        session = Session()
        _dict_id = _get_dict_cache_for_session(get_default_session())
        close_session=True
    else:
        _dict_id = _get_dict_cache_for_session(session)
        close_session=False

    try:
        return _dict_id[text]
    except KeyError:

        try:
            obj = session.query(DictionaryItem).filter_by(text=text).first()
        except:
            if default is None:
                raise
            else:
                return default
        finally:
            if close_session:
                session.close()

        if obj is None:
            if default is None:
                raise
            else:
                return default

        _dict_id[text] = obj.id
        return obj.id


def get_or_create_dictionary_item(session, name):
    """This tries to get the DictionaryItem corresponding to name from
    the database.  If it doesn't exist, it creates a pending
    object. Note that this must be called *while the database is
    locked under the specified session* to prevent duplicate items
    being created"""

    if session not in _dict_obj:
        _dict_obj[session] = {}

    # try to get it from the cache
    obj = _dict_obj[session].get(name, None)

    if obj is not None:
        return obj

    # try to get it from the db
    obj = session.query(DictionaryItem).filter_by(text=name).first()


    if obj is None:
        # try to create it
        try:
            obj = DictionaryItem(name)
            obj = session.merge(obj)
            session.commit()
        except sqlalchemy.exc.IntegrityError:
            session.rollback()
            obj = session.query(DictionaryItem).filter_by(text=name).first()
            if obj is None:
                raise # can't get it from the DB, can't create it from the DB... who knows...

    _dict_obj[session][name] = obj
    return obj

def _get_dict_cache_for_session(session):
    session_dict = _dict_id.get(session, {})
    _dict_id[session] = session_dict
    return session_dict

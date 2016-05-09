# -*- coding: utf-8 -*-
"""
    sphinx.ext.doxytag
    ~~~~~~~~~~~~~~~~~~

    Link to named tags in a Doxygen tagfile.

    :copyright: Copyright 2007-2016 by the Sphinx team, see AUTHORS. Pierre de
        Buyl for this extension, based on todo.py and extlinks.py from the
        Sphinx distribution.
    :license: BSD, see LICENSE for details.
"""

from docutils import nodes, utils
from docutils.parsers.rst import directives

import sphinx
from sphinx.locale import _
from sphinx.util.nodes import split_explicit_title

import os.path
from xml.etree import ElementTree

def load_tagfile(app):
    tagfile = ElementTree.parse(os.path.join(app.srcdir, app.config.doxytag_tagfile))
    env = app.builder.env
    d = {}
    for f in tagfile.findall('compound'):
        if f.get('kind') == 'file':
            for m in f.findall('member'):
                name = m.find('name').text
                anchor_file = m.find('anchorfile').text
                anchor = m.find('anchor').text
                if name not in d:
                    d[name] = anchor_file+'#'+anchor

    env.doxytag_dict = d

def get_doxytag_role(app):
    def doxytag_role(typ, rawtext, text, lineno, inliner, options={}, content=[]):
        env = app.builder.env
        config = app.config
        text = utils.unescape(text)
        full_url = config.doxytag_relative_path + env.doxytag_dict[text]
        pnode = nodes.reference(text, text, internal=False, refuri=full_url)
        return [pnode], []
    return doxytag_role

def setup(app):
    app.add_config_value('doxytag_tagfile', '', 'env')
    app.add_config_value('doxytag_relative_path', '.', 'env')
    app.connect('builder-inited', load_tagfile)
    app.add_role('doxytag', get_doxytag_role(app))

    return {'version': sphinx.__display_version__, 'parallel_read_safe': True}

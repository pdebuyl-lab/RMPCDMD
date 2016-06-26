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

import docutils
from docutils import nodes, utils
from docutils.parsers.rst import directives
from docutils.parsers.rst import Directive

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

def process_line(line):
    line = line[2:].rstrip()
    line = line.lstrip()
    if r'\param' in line:
        line_s = line.split()
        line = '- **'+line_s[1]+'** ' + ' '.join(line_s[2:])
    return line

class DoxyHeaderDirective(Directive):
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = True
    has_content = False

    def run(self):
        env = self.state.document.settings.env
        filename = self.arguments[0]
        data = []
        in_header = False
        in_params = False
        search_program = False
        program_name = ''
        file_path = self.arguments[0]
        with open(file_path, 'r') as p_file:
            for line in p_file:
                if not in_header:
                    if line[:2]=='!>':
                        in_header = True
                    else:
                        continue
                else:
                    if len(line) >= 2:
                        if line[:2] != '!!':
                            search_program = True
                    else:
                        continue
                if r'\param' in line:
                    if not in_params:
                        in_params = True
                        data.append('')
                        data.append('Parameters')
                        data.append('^^^^^^^^^^')
                        data.append('')

                if not search_program:
                    if in_params and not r'\param' in line:
                        line = '!!' + line[8:]
                    data.append(process_line(line))
                else:
                    if 'program' in line:
                        s_line = line.lstrip().split()
                        if len(s_line) >= 2 and s_line[0] == 'program':
                            program_name = s_line[1]
                            break

        t_data = ['``'+program_name+'``']
        t_data += ['-' * (len(program_name)+4)]
        t_data += ['']
        t_data += ['Synopsis']
        t_data += ['^^^^^^^^']
        t_data += ['']
        t_data += ['Source code: :doxytag:`'+program_name+'`']
        t_data += ['']
        data = t_data + data
        self.state_machine.insert_input(data, 'h')
        return []

def setup(app):
    app.add_config_value('doxytag_tagfile', '', 'env')
    app.add_config_value('doxytag_relative_path', '.', 'env')
    app.connect('builder-inited', load_tagfile)
    app.add_role('doxytag', get_doxytag_role(app))
    app.add_directive('doxyheader', DoxyHeaderDirective)

    return {'version': sphinx.__display_version__, 'parallel_read_safe': True}

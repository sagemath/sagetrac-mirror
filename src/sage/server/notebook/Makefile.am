# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk

notebookdir = $(pythondir)/sage/server/notebook
notebook_PYTHON = \
    all.py applet.py avatars.py cell.py colorize.py config.py conf.py css.py \
    docHTMLProcessor.py gnutls_socket_ssl.py guard.py __init__.py interact.py \
    jquery.py js.py keyboards.py mailsender.py nodoctest.py notebook_object.py \
    notebook.py register.py run_notebook.py sage_email.py sagetex.py secureuser.py \
    server_conf.py smtpsend.py template.py tutorial.py twist.py user_conf.py \
    user_db.py user.py worksheet_conf.py worksheet.py

notebook_compressdir = $(notebookdir)/compress
notebook_compress_PYTHON = \
    compress/all.py \
    compress/BaseConvert.py \
    compress/__init__.py \
    compress/JavaScriptCompressor.py \
    compress/nodoctest.py \
    compress/SourceMap.py

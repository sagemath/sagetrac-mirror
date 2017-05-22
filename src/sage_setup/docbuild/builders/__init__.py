def build_other_doc(args):
    document = args[0]
    name = args[1]
    kwds = args[2]
    args = args[3:]
    logger.warning("\nBuilding %s.\n" % document)
    getattr(get_builder(document), name)(*args, **kwds)


class AllBuilder(object):
    """
    A class used to build all of the documentation.
    """
    def __getattr__(self, attr):
        """
        For any attributes not explicitly defined, we just go through
        all of the documents and call their attr.  For example,
        'AllBuilder().json()' will go through all of the documents
        and call the json() method on their builders.
        """
        from functools import partial
        return partial(self._wrapper, attr)

    def _wrapper(self, name, *args, **kwds):
        """
        This is the function which goes through all of the documents
        and does the actual building.
        """
        import time
        start = time.time()
        docs = self.get_all_documents()
        refs = [x for x in docs if x.endswith('reference')]
        others = [x for x in docs if not x.endswith('reference')]

        # Build the reference manual twice to resolve references.  That is,
        # build once with the inventory builder to construct the intersphinx
        # inventory files, and then build the second time for real.  So the
        # first build should be as fast as possible;
        logger.warning("\nBuilding reference manual, first pass.\n")
        for document in refs:
            getattr(get_builder(document), 'inventory')(*args, **kwds)

        logger.warning("Building reference manual, second pass.\n")
        for document in refs:
            getattr(get_builder(document), name)(*args, **kwds)

        # build the other documents in parallel
        L = [(doc, name, kwds) + args for doc in others]
        build_many(build_other_doc, L)
        logger.warning("Elapsed time: %.1f seconds."%(time.time()-start))
        logger.warning("Done building the documentation!")

    def get_all_documents(self):
        """
        Returns a list of all of the documents. A document is a directory within one of
        the language subdirectories of SAGE_DOC_SRC specified by the global LANGUAGES
        variable.

        EXAMPLES::

            sage: from sage_setup.docbuild import AllBuilder
            sage: documents = AllBuilder().get_all_documents()
            sage: 'en/tutorial' in documents
            True
            sage: documents[0] == 'en/reference'
            True
        """
        documents = []
        for lang in LANGUAGES:
            for document in os.listdir(os.path.join(SAGE_DOC_SRC, lang)):
                if (document not in OMIT
                    and os.path.isdir(os.path.join(SAGE_DOC_SRC, lang, document))):
                    documents.append(os.path.join(lang, document))

        # Ensure that the reference guide is compiled first so that links from
        # the other documents to it are correctly resolved.
        if 'en/reference' in documents:
            documents.remove('en/reference')
        documents.insert(0, 'en/reference')

        return documents


if NUM_THREADS > 1:
    def build_many(target, args):
        from multiprocessing import Pool
        pool = Pool(NUM_THREADS, maxtasksperchild=1)
        # map_async handles KeyboardInterrupt correctly. Plain map and
        # apply_async does not, so don't use it.
        x = pool.map_async(target, args, 1)
        try:
            ret = x.get(99999)
            pool.close()
            pool.join()
        except Exception:
            pool.terminate()
            if ABORT_ON_ERROR:
                raise
        return ret
else:
    def build_many(target, args):
        results = []

        for arg in args:
            try:
                results.append(target(arg))
            except Exception:
                if ABORT_ON_ERROR:
                    raise

        return results


def builder_helper(type):
    """
    Returns a function which builds the documentation for
    output type type.
    """
    def f(self, *args, **kwds):
            output_dir = self._output_dir(type)

            options = ALLSPHINXOPTS

            if self.name == 'website':
                # WEBSITESPHINXOPTS is either empty or " -A hide_pdf_links=1 "
                options += WEBSITESPHINXOPTS

            if kwds.get('use_multidoc_inventory', True):
                options += ' -D multidoc_first_pass=0'
            else:
                options += ' -D multidoc_first_pass=1'

            build_command = '-b %s -d %s %s %s %s'%(type, self._doctrees_dir(),
                                                      options, self.dir,
                                                      output_dir)
            logger.debug(build_command)

            # Run Sphinx with Sage's special logger
            sys.argv = ["sphinx-build"] + build_command.split()
            from .sphinxbuild import runsphinx
            try:
                runsphinx()
            except Exception:
                if ABORT_ON_ERROR:
                    raise

            if "/latex" in output_dir:
                logger.warning("LaTeX file written to {}".format(output_dir))
            else:
                logger.warning(
                    "Build finished.  The built documents can be found in {}".
                    format(output_dir))

    f.is_output_format = True
    return f

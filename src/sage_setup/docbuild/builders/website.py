class WebsiteBuilder(DocBuilder):
    def html(self):
        """
        After we've finished building the website index page, we copy
        everything one directory up.  Then we call
        :meth:`create_html_redirects`.
        """
        DocBuilder.html(self)
        html_output_dir = self._output_dir('html')
        for f in os.listdir(html_output_dir):
            src = os.path.join(html_output_dir, f)
            dst = os.path.join(html_output_dir, '..', f)
            if os.path.isdir(src):
                shutil.rmtree(dst, ignore_errors=True)
                shutil.copytree(src, dst)
            else:
                shutil.copy2(src, dst)
        self.create_html_redirects()

    def create_html_redirects(self):
        """
        Writes a number of small HTML files; these are files which
        used to contain the main content of the reference manual
        before before splitting the manual into multiple
        documents. After the split, those files have moved, so in each
        old location, write a file which redirects to the new version.
        (This is so old URLs to pieces of the reference manual still
        open the correct files.)
        """
        # The simple html template which will cause a redirect to the
        # correct file
        html_template = """<html><head>
            <meta HTTP-EQUIV="REFRESH" content="0; url=%s">
            </head><body></body></html>"""

        reference_dir = os.path.abspath(os.path.join(self._output_dir('html'),
                                                     '..', 'reference'))
        reference_builder = ReferenceBuilder('reference')
        refdir = os.path.join(os.environ['SAGE_DOC_SRC'], 'en', 'reference')
        for document in reference_builder.get_all_documents(refdir):
            #path is the directory above reference dir
            path = os.path.abspath(os.path.join(reference_dir, '..'))

            # the name of the subdocument
            document_name = document.split('/')[1]

            # the sage directory within a subdocument, for example
            # local/share/doc/sage/html/en/reference/algebras/sage
            sage_directory = os.path.join(path, document, 'sage')

            # Walk through all of the files in the sage_directory
            for dirpath, dirnames, filenames in os.walk(sage_directory):
                # a string like reference/algebras/sage/algebras
                short_path = dirpath[len(path)+1:]

                # a string like sage/algebras
                shorter_path = os.path.join(*short_path.split(os.sep)[2:])

                #Make the shorter path directory
                try:
                    os.makedirs(os.path.join(reference_dir, shorter_path))
                except OSError:
                    pass

                for filename in filenames:
                    if not filename.endswith('html'):
                        continue

                    # the name of the html file we are going to create
                    redirect_filename = os.path.join(reference_dir, shorter_path, filename)

                    # the number of levels up we need to use in the relative url
                    levels_up = len(shorter_path.split(os.sep))

                    # the relative url that we will redirect to
                    redirect_url = "/".join(['..']*levels_up + [document_name, shorter_path, filename])

                    # write the html file which performs the redirect
                    with open(redirect_filename, 'w') as f:
                        f.write(html_template % redirect_url)


    def clean(self):
        """
        When we clean the output for the website index, we need to
        remove all of the HTML that were placed in the parent
        directory.
        """
        html_output_dir = self._output_dir('html')
        parent_dir = os.path.realpath(os.path.join(html_output_dir, '..'))
        for filename in os.listdir(html_output_dir):
            parent_filename = os.path.join(parent_dir, filename)
            if not os.path.exists(parent_filename):
                continue
            if os.path.isdir(parent_filename):
                shutil.rmtree(parent_filename, ignore_errors=True)
            else:
                os.unlink(parent_filename)

        DocBuilder.clean(self)


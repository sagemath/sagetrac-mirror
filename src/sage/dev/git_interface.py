r"""
Git Interface

This module provides a python interface to Sage's git repository.

AUTHORS:

- David Roe, Julian Rueth, Keshav Kini, Nicolas M. Thiery, Robert Bradshaw:
  initial version

"""
#*****************************************************************************
#       Copyright (C) 2013 David Roe <roed.math@gmail.com>
#                          Julian Rueth <julian.rueth@fsfe.org>
#                          Keshav Kini <keshav.kini@gmail.com>
#                          Nicolas M. Thiery <Nicolas.Thiery@u-psud.fr>
#                          Robert Bradshaw <robertwb@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os

from sage.env import SAGE_DOT_GIT, SAGE_REPO_AUTHENTICATED, SAGE_ROOT

from git_error import GitError, DetachedHeadError

class GitProxy(object):
    r"""
    A proxy object to wrap actual calls to git.

    EXAMPLES::

        sage: from sage.dev.git_interface import GitProxy
        sage: from sage.dev.test.config import DoctestConfig
        sage: from sage.dev.test.user_interface import DoctestUserInterface
        sage: config = DoctestConfig()
        sage: GitProxy(config['git'], DoctestUserInterface(config['UI']))
        <sage.dev.git_interface.GitProxy object at 0x...>

    """
    def __init__(self, config, UI, upload_ssh_key=None):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.git_interface import GitProxy
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: type(GitProxy(config['git'], DoctestUserInterface(config['UI'])))
            <class 'sage.dev.git_interface.GitProxy'>

        """
        self._config = config
        self._UI = UI
        self.__upload_ssh_key = upload_ssh_key

        self._src = self._config.get('src', SAGE_ROOT)
        self._dot_git = self._config.get('dot_git', SAGE_DOT_GIT)
        self._gitcmd = self._config.get('gitcmd', 'git')
        self._repository = self._config.get('repository', SAGE_REPO_AUTHENTICATED)

        if not os.path.isabs(self._src):
            raise ValueError("`%s` is not an absolute path."%self._src)
        if not os.path.exists(self._src):
            raise ValueError("`%s` does not point to an existing directory."%self._src)
        if not os.path.isabs(self._dot_git):
            raise ValueError("`%s` is not an absolute path."%self._dot_git)
        if not os.path.exists(self._dot_git):
            raise ValueError("`%s` does not point to an existing directory."%self._dot_git)

    def _run_git(self, cmd, args, kwds, **ckwds):
        r"""
        Common implementation for :meth:`_execute`, :meth:`_execute_silent`,
        :meth:`_execute_supersilent`, and :meth:`_read_output`

        INPUT:

        - ``cmd`` - git command run

        - ``args`` - extra arguments for git

        - ``kwds`` - extra keywords for git

        - ``ckwds`` - Popen like keywords but with the following changes

          - ``stdout`` - if set to ``False`` will supress stdout

          - ``stderr`` - if set to ``False`` will supress stderr

        .. WARNING::

            This method does not raise an exception if the git call returns a
            non-zero exit code.

        EXAMPLES::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))
            sage: os.chdir(config['git']['src'])

            sage: git._run_git('status', (), {})
            # On branch master
            #
            # Initial commit
            #
            nothing to commit (create/copy files and use "git add" to track)
            (0, None, None, 'git -c user.email=doc@test.test -c user.name=doctest status')
            sage: git._run_git('status', (), {}, stdout=False)
            (0, None, None, 'git -c user.email=doc@test.test -c user.name=doctest status')

        TESTS:

        Check that we refuse to touch the live source code in doctests::

            sage: dev.git.status()
            Traceback (most recent call last):
            ...
            AssertionError: possible attempt to work with the live repository/directory in a doctest

        """
        import sage.doctest
        import os
        assert not sage.doctest.DOCTEST_MODE or (self._dot_git != SAGE_DOT_GIT and self._repository != SAGE_REPO_AUTHENTICATED and os.path.abspath(self._src).startswith(self._src)), "possible attempt to work with the live repository/directory in a doctest"

        # not sure which commands could possibly create a commit object with
        # unless there are some crazy flags set - these commands should be safe
        if cmd not in [ "config", "diff", "grep", "log", "ls_remote", "remote", "reset", "show", "show_ref", "status", "symbolic_ref" ]:
            self._check_user_email()
        # which commits access the remote repository?
        # again, these should be safe
        self._upload_ssh_key()

        s = [self._gitcmd, "--git-dir=%s"%self._dot_git, "--work-tree=%s"%self._src, cmd]
        if 'user.name' in self._config:
            s.insert(3, '-c')
            s.insert(4, 'user.name='+self._config['user.name'])
        if 'user.email' in self._config:
            s.insert(3, '-c')
            s.insert(4, 'user.email='+self._config['user.email'])

        env = ckwds.setdefault('env', dict(os.environ))
        env.update(kwds.pop('env', {}))

        for k, v in kwds.iteritems():
            if len(k) == 1:
                k = '-' + k
            else:
                k = '--' + k.replace('_', '-')
            if v is True:
                s.append(k)
            elif v is not False:
                s.extend((k, v))
        if args:
            s.extend(a for a in args if a is not None)

        s = [str(arg) for arg in s]

        complete_cmd = " ".join([arg for i,arg in enumerate(s) if i not in (1,2)]) # drop --git-dir, --work-tree from debug output
        self._UI.debug("[git] %s"%complete_cmd)

        if ckwds.get('dryrun', False):
            return s

        import subprocess
        drop_stdout = ckwds.get('stdout') is False
        read_stdout = ckwds.get('stdout') is str
        drop_stderr = ckwds.get('stderr') is False
        read_stderr = ckwds.get('stderr') is str

        if drop_stdout or read_stdout:
            ckwds['stdout'] = subprocess.PIPE
        if drop_stderr or read_stderr:
            ckwds['stderr'] = subprocess.PIPE

        process = subprocess.Popen(s, **ckwds)
        stdout, stderr = process.communicate()
        retcode = process.poll()

        # recover stdout and stderr for debugging on non-zero exit code
        if retcode:
            if drop_stdout or read_stdout:
                pass
            else:
                stdout = None

            if drop_stderr or read_stderr:
                pass
            else:
                stderr = None
        else:
            if not read_stdout:
                stdout = None
            if not read_stderr:
                stderr = None

        return retcode, stdout, stderr, complete_cmd

    def _execute(self, cmd, *args, **kwds):
        r"""
        Run git.

        Raises an exception if git has non-zero exit code.

        INPUT:

        - ``cmd`` - git command run

        - ``args`` - extra arguments for git

        - ``kwds`` - extra keywords for git

        EXAMPLES::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))
            sage: os.chdir(config['git']['src'])

            sage: git._execute('status')
            # On branch master
            #
            # Initial commit
            #
            nothing to commit (create/copy files and use "git add" to track)
            sage: git._execute('status',foo=True) # --foo is not a valid parameter
            Traceback (most recent call last):
            ...
            GitError: git returned with non-zero exit code (129) for `git -c user.email=doc@test.test -c user.name=doctest status --foo`.

        """
        exit_code, stdout, stderr, cmd = self._run_git(cmd, args, kwds)
        if exit_code:
            raise GitError(exit_code, cmd, stdout, stderr)

    def _execute_silent(self, cmd, *args, **kwds):
        r"""
        Run git and supress its output to stdout.

        Same input as :meth:`execute`.

        Raises an error if git returns a non-zero exit code.

        EXAMPLES::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))
            sage: os.chdir(config['git']['src'])

            sage: git._execute_silent('status')
            sage: git._execute_silent('status',foo=True) # --foo is not a valid parameter
            Traceback (most recent call last):
            ...
            GitError: git returned with non-zero exit code (129) for `git -c user.email=doc@test.test -c user.name=doctest status --foo`.

        """
        exit_code, stdout, stderr, cmd = self._run_git(cmd, args, kwds, stdout=False)
        if exit_code:
            raise GitError(exit_code, cmd, stdout, stderr)

    def _execute_supersilent(self, cmd, *args, **kwds):
        r"""
        Run git and supress its output to stdout and stderr.

        Same input as :meth:`execute`.

        Raises an error if git returns a non-zero exit code.

        EXAMPLES::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))
            sage: os.chdir(config['git']['src'])

            sage: git._execute_supersilent('status')
            sage: git._execute_supersilent('status',foo=True) # --foo is not a valid parameter
            Traceback (most recent call last):
            ...
            GitError: git returned with non-zero exit code (129) for `git -c user.email=doc@test.test -c user.name=doctest status --foo`.
            ...

        """
        exit_code, stdout, stderr, cmd = self._run_git(cmd, args, kwds, stdout=False, stderr=False)
        if exit_code:
            raise GitError(exit_code, cmd, stdout, stderr)

    def _read_output(self, cmd, *args, **kwds):
        r"""
        Run git and return its output to stdout.

        Same input as :meth:`execute`.

        Raises an error if git returns a non-zero exit code.

        EXAMPLES::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))
            sage: os.chdir(config['git']['src'])

            sage: git._read_output('status')
            '# On branch master\n#\n# Initial commit\n#\nnothing to commit (create/copy files and use "git add" to track)\n'
            sage: git._read_output('status',foo=True) # --foo is not a valid parameter
            Traceback (most recent call last):
            ...
            GitError: git returned with non-zero exit code (129) for `git -c user.email=doc@test.test -c user.name=doctest status --foo`.
            ...

        """
        exit_code, stdout, stderr, cmd = self._run_git(cmd, args, kwds, stdout=str, stderr=False)
        if exit_code:
            raise GitError(exit_code, cmd, stdout, stderr)
        return stdout

    def _check_user_email(self):
        r"""
        Make sure that a real name and an email are set for git. These will
        show up next to any commit that user creates.

        TESTS::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: del config['git']['user.name']
            sage: del config['git']['user.email']
            sage: UI = DoctestUserInterface(config["UI"])
            sage: git = GitInterface(config["git"], UI)
            sage: os.chdir(config['git']['src'])
            sage: UI.append("Doc Test")
            sage: UI.append("doc@test")
            sage: git._check_user_email() # random output, because it depends on whether the user has set user.name and user.email in its .gitconfig

        """
        if self._config.get('user_email_set', False):
            return

        if 'user.name' not in self._config:
            try:
                self._execute_supersilent("config","user.name")
            except GitError as e:
                if e.exit_code == 1:
                    name = self._UI.get_input("No real name has been set for git. This name shows up as the author for any commits you contribute to sage. Your real name:")
                    self._execute("config","user.name",name,local=True,add=True)
                    self._UI.info("Your real name has been saved.")
                else:
                    raise

        if 'user.email' not in self._config:
            try:
                self._execute_supersilent("config", "user.email")
            except GitError as e:
                if e.exit_code == 1:
                    email = self._UI._get_input("No email address has been set for git. This email shows up as the author for any commits you contribute to sage. Your email address:")
                    self._execute("config","user.email",email,local=True,add=True)
                    self._UI.info("Your email has been saved.")
                else:
                    raise

        self._config['user_email_set'] = "True"

    def _upload_ssh_key(self):
        r"""
        Make sure that the public ssh key has been uploaded to the trac server.

        TESTS::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: del config['git']['user.name']
            sage: del config['git']['user.email']
            sage: UI = DoctestUserInterface(config["UI"])
            sage: def upload_ssh_key(): print "Uploading ssh key."
            sage: git = GitInterface(config["git"], UI, upload_ssh_key=upload_ssh_key)
            sage: os.chdir(config['git']['src'])
            sage: from sage.dev.git_error import GitError
            sage: try:
            ....:   git.ls_remote()
            ....: except GitError: pass
            Uploading ssh key.
            sage: try:
            ....:   git.ls_remote()
            ....: except GitError: pass

        """
        if self._config.get('ssh_key_set', False):
            return

        if self.__upload_ssh_key is not None:
            from user_interface_error import OperationCancelledError
            try:
                self.__upload_ssh_key()
            except OperationCancelledError:
                pass # do not bother the user again, probably the key has been uploaded manually already

            self._config['ssh_key_set'] = "True"

class ReadStdoutGitProxy(GitProxy):
    r"""
    A proxy object to wrap calls to git.

    Calls to git return the stdout of git and raise an error on a non-zero exit
    code. Output to stderr is supressed.

    EXAMPLES::

        sage: dev.git.status() # not tested

    """
    __call__ = GitProxy._read_output

class SilentGitProxy(GitProxy):
    r"""
    A proxy object to wrap calls to git.

    Calls to git do not show any output to stdout and raise an error on a
    non-zero exit code. Output to stderr is printed.

    EXAMPLES::

        sage: dev.git.silent.status() # not tested

    """
    __call__ = GitProxy._execute_silent

class EchoGitProxy(GitProxy):
    r"""
    A proxy object to wrap calls to git.

    Calls to git show output to stdout and stderr as if the command was
    executed directly and raise an error on a non-zero exit code.

    EXAMPLES::

        sage: dev.git.echo.status() # not tested

    """
    __call__ = GitProxy._execute

class SuperSilentGitProxy(GitProxy):
    r"""
    A proxy object to wrap calls to git.

    Calls to git do not show any output to stderr or stdout and raise an error
    on a non-zero exit code.

    EXAMPLES::

        sage: dev.git.super_silent.status() # not tested

    """
    __call__ = GitProxy._execute_supersilent

class GitInterface(ReadStdoutGitProxy):
    r"""
    A wrapper around the ``git`` command line tool.

    Most methods of this class correspond to actual git commands. Some add
    functionality which is not directly available in git. However, all of the
    methods should be non-interactive. If interaction is required the method
    should live in :class:`saged.dev.sagedev.SageDev`.

    EXAMPLES::

        sage: from sage.dev.config import Config
        sage: from sage.dev.user_interface import UserInterface
        sage: from sage.dev.git_interface import GitInterface
        sage: config = Config()
        sage: GitInterface(config['git'], UserInterface(config['UI']))
        GitInterface()

    """
    def __init__(self, config, UI, upload_ssh_key=None):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.config import Config
            sage: from sage.dev.user_interface import UserInterface
            sage: from sage.dev.git_interface import GitInterface
            sage: config = Config()
            sage: type(GitInterface(config['git'], UserInterface(config['UI'])))
            <class 'sage.dev.git_interface.GitInterface'>

        """
        ReadStdoutGitProxy.__init__(self, config, UI, upload_ssh_key=upload_ssh_key)

        self.silent = SilentGitProxy(config, UI, upload_ssh_key=upload_ssh_key)
        self.super_silent = SuperSilentGitProxy(config, UI, upload_ssh_key=upload_ssh_key)
        self.echo = EchoGitProxy(config, UI, upload_ssh_key=upload_ssh_key)

    def __repr__(self):
        r"""
        Return a printable representation of this object.

        TESTS::

            sage: from sage.dev.config import Config
            sage: from sage.dev.user_interface import UserInterface
            sage: from sage.dev.git_interface import GitInterface
            sage: config = Config()
            sage: repr(GitInterface(config['git'], UserInterface(config['UI'])))
            'GitInterface()'

        """
        return "GitInterface()"

    def get_state(self):
        r"""
        Get the current state of merge/rebase/am/etc operations.

        OUTPUT:

        A tuple of strings which consists of any of the following:
        ``'rebase'``, ``'am'``, ``'rebase-i'``, ``'rebase-m'``, ``'merge'``,
        ``'bisect'``, ``'cherry-seq'``, ``'cherry'``.

        EXAMPLES:

        Create a :class:`GitInterface` for doctesting::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        Create two conflicting branches::

            sage: os.chdir(config['git']['src'])
            sage: with open("file","w") as f: f.write("version 0")
            sage: git.silent.add("file")
            sage: git.silent.commit("-m","initial commit")
            sage: git.super_silent.checkout("-b","branch1")
            sage: with open("file","w") as f: f.write("version 1")
            sage: git.silent.commit("-am","second commit")
            sage: git.super_silent.checkout("master")
            sage: git.super_silent.checkout("-b","branch2")
            sage: with open("file","w") as f: f.write("version 2")
            sage: git.silent.commit("-am","conflicting commit")

        A ``merge`` state::

            sage: git.super_silent.checkout("branch1")
            sage: git.silent.merge('branch2')
            Traceback (most recent call last):
            ...
            GitError: git returned with non-zero exit code (1) for `git -c user.email=doc@test.test -c user.name=doctest merge branch2`.
            ...
            sage: git.get_state()
            ('merge',)
            sage: git.silent.merge(abort=True)
            sage: git.get_state()
            ()

        A ``rebase`` state::

            sage: git._execute_supersilent('rebase', 'branch2')
            Traceback (most recent call last):
            ...
            GitError: git returned with non-zero exit code (1) for `git -c user.email=doc@test.test -c user.name=doctest rebase branch2`.
            ...
            sage: git.get_state()
            ('rebase',)
            sage: git.super_silent.rebase(abort=True)
            sage: git.get_state()
            ()

        A merge within an interactive rebase::

            sage: git.super_silent.rebase('HEAD^', interactive=True, env={'GIT_SEQUENCE_EDITOR':'sed -i s+pick+edit+'})
            sage: git.get_state()
            ('rebase-i',)
            sage: git.silent.merge('branch2')
            Traceback (most recent call last):
            ...
            GitError: git returned with non-zero exit code (1) for `git -c user.email=doc@test.test -c user.name=doctest merge branch2`.
            ...
            sage: git.get_state()
            ('merge', 'rebase-i')
            sage: git.super_silent.rebase(abort=True)
            sage: git.get_state()
            ()

        """
        # logic based on zsh's git backend for vcs_info
        opj = os.path.join
        p = lambda x: opj(self._dot_git, x)
        ret = []
        for d in map(p,("rebase-apply", "rebase", opj("..",".dotest"))):
            if os.path.isdir(d):
                if os.path.isfile(opj(d, 'rebasing')) and 'rebase' not in ret:
                    ret.append('rebase')
                if os.path.isfile(opj(d, 'applying')) and 'am' not in ret:
                    ret.append('am')
        for f in map(p, (opj('rebase-merge', 'interactive'),
                         opj('.dotest-merge', 'interactive'))):
            if os.path.isfile(f):
                ret.append('rebase-i')
                break
        else:
            for d in map(p, ('rebase-merge', '.dotest-merge')):
                if os.path.isdir(d):
                    ret.append('rebase-m')
                    break
        if os.path.isfile(p('MERGE_HEAD')):
            ret.append('merge')
        if os.path.isfile(p('BISECT_LOG')):
            ret.append('bisect')
        if os.path.isfile(p('CHERRY_PICK_HEAD')):
            if os.path.isdir(p('sequencer')):
                ret.append('cherry-seq')
            else:
                ret.append('cherry')
        # return in reverse order so reset operations are correctly ordered
        return tuple(reversed(ret))

    def reset_to_clean_state(self):
        r"""
        Get out of a merge/am/rebase/etc state.

        EXAMPLES:

        Create a :class:`GitInterface` for doctesting::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        Create two conflicting branches::

            sage: os.chdir(config['git']['src'])
            sage: with open("file","w") as f: f.write("version 0")
            sage: git.silent.add("file")
            sage: git.silent.commit("-m","initial commit")
            sage: git.super_silent.checkout("-b","branch1")
            sage: with open("file","w") as f: f.write("version 1")
            sage: git.silent.commit("-am","second commit")
            sage: git.super_silent.checkout("master")
            sage: git.super_silent.checkout("-b","branch2")
            sage: with open("file","w") as f: f.write("version 2")
            sage: git.silent.commit("-am","conflicting commit")

        A merge within an interactive rebase::

            sage: git.super_silent.checkout("branch1")
            sage: git.super_silent.rebase('HEAD^', interactive=True, env={'GIT_SEQUENCE_EDITOR':'sed -i s+pick+edit+'})
            sage: git.get_state()
            ('rebase-i',)
            sage: git.silent.merge('branch2')
            Traceback (most recent call last):
            ...
            GitError: git returned with non-zero exit code (1) for `git -c user.email=doc@test.test -c user.name=doctest merge branch2`.
            ...
            sage: git.get_state()
            ('merge', 'rebase-i')

        Get out of this state::

            sage: git.reset_to_clean_state()
            sage: git.get_state()
            ()

        """
        states = self.get_state()
        if not states:
            return

        state = states[0]
        if state.startswith('rebase'):
            self.silent.rebase(abort=True)
        elif state == 'am':
            self.silent.am(abort=True)
        elif state == 'merge':
            self.silent.merge(abort=True)
        elif state == 'bisect':
            self.silent.bisect(reset=True)
        elif state.startswith('cherry'):
            self.silent.cherry_pick(abort=True)
        else:
            raise RuntimeError("'%s' is not a valid state"%state)

        return self.reset_to_clean_state()

    def reset_to_clean_working_directory(self, remove_untracked_files=False, remove_untracked_directories=False, remove_ignored=False):
        r"""
        Reset any changes made to the working directory.

        INPUT:

        - ``remove_untracked_files`` -- a boolean (default: ``False``), whether
          to remove files which are not tracked by git

        - ``remove_untracked_directories`` -- a boolean (default: ``False``),
          whether to remove directories which are not tracked by git

        - ``remove_ignored`` -- a boolean (default: ``False``), whether to
          remove files directories which are ignored by git

        EXAMPLES:

        Create a :class:`GitInterface` for doctesting::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        Set up some files/directories::

            sage: os.chdir(config['git']['src'])
            sage: open('tracked','w').close()
            sage: git.silent.add('tracked')
            sage: with open('.gitignore','w') as f: f.write('ignored\nignored_dir')
            sage: git.silent.add('.gitignore')
            sage: git.silent.commit('-m', 'initial commit')

            sage: os.mkdir('untracked_dir')
            sage: open('untracked_dir/untracked','w').close()
            sage: open('untracked','w').close()
            sage: open('ignored','w').close()
            sage: os.mkdir('ignored_dir')
            sage: open('ignored_dir/untracked','w').close()
            sage: with open('tracked','w') as f: f.write('version 0')
            sage: git.echo.status()
            # On branch master
            # Changes not staged for commit:
            #   (use "git add <file>..." to update what will be committed)
            #   (use "git checkout -- <file>..." to discard changes in working directory)
            #
            #   modified:   tracked
            #
            # Untracked files:
            #   (use "git add <file>..." to include in what will be committed)
            #
            #   untracked
            #   untracked_dir/
            no changes added to commit (use "git add" and/or "git commit -a")

        Some invalid combinations of flags::

            sage: git.reset_to_clean_working_directory(remove_untracked_files = False, remove_untracked_directories = True)
            Traceback (most recent call last):
            ...
            ValueError: remove_untracked_directories only valid if remove_untracked_files is set
            sage: git.reset_to_clean_working_directory(remove_untracked_files = False, remove_ignored = True)
            Traceback (most recent call last):
            ...
            ValueError: remove_ignored only valid if remove_untracked_files is set

        Per default only the tracked modified files are reset to a clean
        state::

            sage: git.reset_to_clean_working_directory()
            sage: git.echo.status()
            # On branch master
            # Untracked files:
            #   (use "git add <file>..." to include in what will be committed)
            #
            #   untracked
            #   untracked_dir/
            nothing added to commit but untracked files present (use "git add" to track)

        Untracked items can be removed by setting the parameters::

            sage: git.reset_to_clean_working_directory(remove_untracked_files=True)
            Removing untracked
            Not removing untracked_dir/
            sage: git.reset_to_clean_working_directory(remove_untracked_files=True, remove_untracked_directories=True)
            Removing untracked_dir/
            sage: git.reset_to_clean_working_directory(remove_untracked_files=True, remove_ignored=True)
            Removing ignored
            Not removing ignored_dir/
            sage: git.reset_to_clean_working_directory(remove_untracked_files=True, remove_untracked_directories=True, remove_ignored=True)
            Removing ignored_dir/

        """
        if remove_untracked_directories and not remove_untracked_files:
            raise ValueError("remove_untracked_directories only valid if remove_untracked_files is set")
        if remove_ignored and not remove_untracked_files:
            raise ValueError("remove_ignored only valid if remove_untracked_files is set")

        self.silent.reset(hard=True)

        if remove_untracked_files:
            switches = ['-f']
            if remove_untracked_directories: switches.append("-d")
            if remove_ignored: switches.append("-x")
            self.echo.clean(*switches)

    def is_child_of(self, a, b):
        r"""
        Return whether ``a`` is a child of ``b``.

        EXAMPLES:

        Create a :class:`GitInterface` for doctesting::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        Create two conflicting branches::

            sage: os.chdir(config['git']['src'])
            sage: with open("file","w") as f: f.write("version 0")
            sage: git.silent.add("file")
            sage: git.silent.commit("-m","initial commit")
            sage: git.super_silent.checkout("-b","branch1")
            sage: with open("file","w") as f: f.write("version 1")
            sage: git.silent.commit("-am","second commit")
            sage: git.super_silent.checkout("master")
            sage: git.super_silent.checkout("-b","branch2")
            sage: with open("file","w") as f: f.write("version 2")
            sage: git.silent.commit("-am","conflicting commit")

            sage: git.is_child_of('master', 'branch2')
            False
            sage: git.is_child_of('branch2', 'master')
            True
            sage: git.is_child_of('branch1', 'branch2')
            False
            sage: git.is_child_of('master', 'master')
            True

        """
        return self.is_ancestor_of(b, a)

    def is_ancestor_of(self, a, b):
        r"""
        Return whether ``a`` is an ancestor of ``b``.

        EXAMPLES:

        Create a :class:`GitInterface` for doctesting::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        Create two conflicting branches::

            sage: os.chdir(config['git']['src'])
            sage: with open("file","w") as f: f.write("version 0")
            sage: git.silent.add("file")
            sage: git.silent.commit("-m","initial commit")
            sage: git.super_silent.checkout("-b","branch1")
            sage: with open("file","w") as f: f.write("version 1")
            sage: git.silent.commit("-am","second commit")
            sage: git.super_silent.checkout("master")
            sage: git.super_silent.checkout("-b","branch2")
            sage: with open("file","w") as f: f.write("version 2")
            sage: git.silent.commit("-am","conflicting commit")

            sage: git.is_ancestor_of('master', 'branch2')
            True
            sage: git.is_ancestor_of('branch2', 'master')
            False
            sage: git.is_ancestor_of('branch1', 'branch2')
            False
            sage: git.is_ancestor_of('master', 'master')
            True

        """
        return self.merge_base(a, b) == self.rev_parse(a)

    def has_uncommitted_changes(self):
        r"""
        Return whether there are uncommitted changes, i.e., whether there are
        modified files which are tracked by git.

        EXAMPLES:

        Create a :class:`GitInterface` for doctesting::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        An untracked file does not count towards uncommited changes::

            sage: os.chdir(config['git']['src'])
            sage: open('untracked','w').close()
            sage: git.has_uncommitted_changes()
            False

        Once added to the index it does::

            sage: git.silent.add('untracked')
            sage: git.has_uncommitted_changes()
            True
            sage: git.silent.commit('-m', 'tracking untracked')
            sage: git.has_uncommitted_changes()
            False
            sage: with open('untracked','w') as f: f.write('version 0')
            sage: git.has_uncommitted_changes()
            True

        """
        return bool([line for line in self.status(porcelain=True).splitlines() if not line.startswith('?')])

    def untracked_files(self):
        r"""
        Return a list of file names for files that are not tracked by git and
        not ignored.

        EXAMPLES:

        Create a :class:`GitInterface` for doctesting::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        An untracked file::

            sage: os.chdir(config['git']['src'])
            sage: git.untracked_files()
            []
            sage: open('untracked','w').close()
            sage: git.untracked_files()
            ['untracked']

         Directories are not displayed here::

            sage: os.mkdir('untracked_dir')
            sage: git.untracked_files()
            ['untracked']
            sage: open('untracked_dir/untracked','w').close()
            sage: git.untracked_files()
            ['untracked', 'untracked_dir/untracked']

        """
        import os
        old_cwd = os.getcwd()
        if 'src' in self._config:
            os.chdir(self._config['src'])
        else:
            from sage.env import SAGE_ROOT
            os.chdir(SAGE_ROOT)
        try:
            fnames = self.ls_files(other=True, exclude_standard=True).splitlines()
            fnames = [ os.path.abspath(fname) for fname in fnames ]
            return [ os.path.relpath(fname, old_cwd) for fname in fnames ]
        finally:
            os.chdir(old_cwd)

    def local_branches(self):
        r"""
        Return a list of local branches sorted by last commit time.

        EXAMPLES::

        Create a :class:`GitInterface` for doctesting::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        Create some branches::

            sage: os.chdir(config['git']['src'])
            sage: git.silent.commit('-m','initial commit','--allow-empty')
            sage: git.silent.branch('branch1')
            sage: git.silent.branch('branch2')

        Use this repository as a remote repository::

            sage: config2 = DoctestConfig()
            sage: git2 = GitInterface(config2["git"], DoctestUserInterface(config["UI"]))
            sage: os.chdir(config2['git']['src'])
            sage: git2.silent.commit('-m','initial commit','--allow-empty')
            sage: git2.silent.remote('add', 'git', config['git']['src'])
            sage: git2.super_silent.fetch('git')
            sage: git2.super_silent.checkout("branch1")
            sage: git2.echo.branch("-a")
            * branch1
              master
              remotes/git/branch1
              remotes/git/branch2
              remotes/git/master

            sage: git2.local_branches()
            ['branch1', 'master']
            sage: os.chdir(config['git']['src'])
            sage: git.local_branches()
            ['branch1', 'branch2', 'master']

        """
        result = self.for_each_ref('refs/heads/', sort='-committerdate', format="%(refname)").splitlines()
        return sorted([head[11:] for head in result])

    def current_branch(self):
        r"""
        Return the current branch

        EXAMPLES:

        Create a :class:`GitInterface` for doctesting::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        Create some branches::

            sage: os.chdir(config['git']['src'])
            sage: git.silent.commit('-m','initial commit','--allow-empty')
            sage: git.silent.commit('-m','second commit','--allow-empty')
            sage: git.silent.branch('branch1')
            sage: git.silent.branch('branch2')

            sage: git.current_branch()
            'master'
            sage: git.super_silent.checkout('branch1')
            sage: git.current_branch()
            'branch1'

        If ``HEAD`` is detached::

            sage: git.super_silent.checkout('master~')
            sage: git.current_branch()
            Traceback (most recent call last):
            ...
            DetachedHeadError: unexpectedly, git is in a detached HEAD state

        """
        try:
            return self.symbolic_ref('HEAD', short=True, quiet=True).strip()
        except GitError as e:
            if e.exit_code == 1:
               raise DetachedHeadError()
            raise

    def commit_for_branch(self, branch):
        r"""
        Return the commit id of the local ``branch``, or ``None`` if the branch
        does not exist

        EXAMPLES:

        Create a :class:`GitInterface` for doctesting::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        Create some branches::

            sage: os.chdir(config['git']['src'])
            sage: git.silent.commit('-m','initial commit','--allow-empty')
            sage: git.silent.branch('branch1')
            sage: git.silent.branch('branch2')

        Check existence of branches::

            sage: git.commit_for_branch('branch1') # random output
            '087e1fdd0fe6f4c596f5db22bc54567b032f5d2b'
            sage: git.commit_for_branch('branch2') is not None
            True
            sage: git.commit_for_branch('branch3') is not None
            False

        """
        return self.commit_for_ref("refs/heads/%s"%branch)

    def commit_for_ref(self, ref):
        r"""
        Return the commit id of the ``ref``, or ``None`` if the ``ref`` does
        not exist.

        EXAMPLES:

        Create a :class:`GitInterface` for doctesting::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        Create some branches::

            sage: os.chdir(config['git']['src'])
            sage: git.silent.commit('-m','initial commit','--allow-empty')
            sage: git.silent.branch('branch1')
            sage: git.silent.branch('branch2')

        Check existence of branches::

            sage: git.commit_for_ref('refs/heads/branch1') # random output
            '087e1fdd0fe6f4c596f5db22bc54567b032f5d2b'
            sage: git.commit_for_ref('refs/heads/branch2') is not None
            True
            sage: git.commit_for_ref('refs/heads/branch3') is not None
            False

        """
        try:
            return self.rev_parse(ref, verify=True).strip()
        except GitError:
            return None

    def rename_branch(self, oldname, newname):
        r"""
        Rename ``oldname`` to ``newname``.

        EXAMPLES:

        Create a :class:`GitInterface` for doctesting::

            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        Create some branches::

            sage: os.chdir(config['git']['src'])
            sage: git.silent.commit('-m','initial commit','--allow-empty')
            sage: git.silent.branch('branch1')
            sage: git.silent.branch('branch2')

        Rename some branches::

            sage: git.rename_branch('branch1', 'branch3')
            sage: git.rename_branch('branch2', 'branch3')
            Traceback (most recent call last):
            ...
            GitError: git returned with non-zero exit code (128) for `git -c user.email=doc@test.test -c user.name=doctest branch --move branch2 branch3`.
            output to stderr: fatal: A branch named 'branch3' already exists.

        """
        self.branch(oldname, newname, move=True)

for git_cmd_ in (
        "add",
        "am",
        "apply",
        "bisect",
        "branch",
        "config",
        "checkout",
        "cherry_pick",
        "clean",
        "clone",
        "commit",
        "diff",
        "fetch",
        "for_each_ref",
        "format_patch",
        "grep",
        "init",
        "log",
        "ls_files",
        "ls_remote",
        "merge",
        "merge_base",
        "mv",
        "pull",
        "push",
        "rebase",
        "remote",
        "reset",
        "rev_list",
        "rev_parse",
        "rm",
        "show",
        "show_ref",
        "stash",
        "status",
        "symbolic_ref",
        "tag"
        ):
    def create_wrapper(git_cmd__):
        r"""
        Create a wrapper for ``git_cmd__``.

        EXAMPLES::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))
            sage: os.chdir(config['git']['src'])
            sage: git.echo.status()
            # On branch master
            #
            # Initial commit
            #
            nothing to commit (create/copy files and use "git add" to track)

        """
        git_cmd = git_cmd__.replace("_","-")
        def meth(self, *args, **kwds):
            return self(git_cmd, *args, **kwds)
        meth.__doc__ = r"""
        Call `git {0}`.

        OUTPUT:

        See the docstring of ``__call__`` for more information.

        EXAMPLES:

            sage: dev.git.{1}() # not tested

        """.format(git_cmd, git_cmd__)
        return meth
    setattr(GitProxy, git_cmd_, create_wrapper(git_cmd_))


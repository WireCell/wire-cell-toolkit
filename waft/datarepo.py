# Hard-coded list of WCT release versions for which THIS release
# should provide historical files.
#
# Guidance for developers and release managers:
#
# - This list should be ammended for each new release.
# - At very least, new release should be added.
# - Include only tagged releases (not "dirty" versions).
# - Pruning of the list is fine.
TEST_DATA_VERSIONS='0.20.0,0.21.0,0.22.0,0.23.0,0.24.1'

import io
from waflib.Logs import debug, info, error, warn
from waflib.Task import Task
try:
    from urllib import request
except ImportError:
    from urllib import urlopen
else:
    urlopen = request.urlopen
import tarfile

import ssl
# www.phy.bnl.gov cert is in a format that is legit but not correctly parsed by python
ssl_context = ssl._create_unverified_context()

def options(opt):
    # fixme: need to provide "official" URL.  For now:
    # https://www.phy.bnl.gov/~bviren/tmp/wcttest/data_repo
    opt.add_option('--test-data-url', type=str, default="https://www.phy.bnl.gov/~bviren/tmp/wcttest/data_repo",
                   help="URL for test data repository [default: see source]")
    opt.add_option('--test-data-versions', type=str, default=TEST_DATA_VERSIONS,
                   help="comma-separated list of tags for which historical test data files are downloaded [default: see source]")


def configure(cfg):
    
    # fixme: MUST load datarepo after wcb_unit_test
    # fixme: refactor datarepo, smplpkgs and wcb_unit_test

    versions = cfg.options.test_data_versions
    versions = versions.split(",")
    cfg.env.TEST_DATA_VERSIONS = versions

    url = cfg.options.test_data_url
    if not url:
        info("datarepo: test data will not be available")
        return
    cfg.env.TEST_DATA_URL = url
    cfg.env.TEST_DATA_DIR = "tests"


def download_and_untar(task):
    '''
    Rule function to download URL as a tar file and unpack under path.
    '''
    unode = task.inputs[0]
    url = unode.read()
    path = unode.parent().abspath()

    debug(f'datarepo: getting {url}')
    web = urlopen(url, context=ssl_context)
    if web.getcode() != 200:
        return web.getcode()
    # tarfile needs seekable file object
    buf = io.BytesIO(web.read())
    buf.seek(0)

    try:
        tf = tarfile.open(fileobj=buf)
        tf.extractall(path)
    except tarfile.TarError as err:
        error(f"datarepo: {err}")
        return -1

class seed_datarepo_urls(Task):
    '''
    Make .url files
    '''

    vars = ['TEST_DATA_VERSIONS', 'TEST_DATA_URL']

    def run(self):
        url = self.env.TEST_DATA_URL
        if not url:
            return
        versions = self.env.TEST_DATA_VERSIONS
        if not versions:
            return
        if isinstance(versions, str):
            versions = versions.split(",")

        bld = self.generator.bld
        tdir = bld.bldnode.make_node(self.env.TEST_DATA_DIR)
        tdir.mkdir()

        for ver in versions:
            uf = tdir.make_node(f'history-{ver}.url')
            uf.write(f'{url}/history-{ver}.tar')

class datarepo_provision(Task):
    vars = ['datarepo_url', 'datarepo_category', 'datarepo_version']
    def run(self):
        
        url = self.env.datarepo_url + "/" + self.env.datarepo_category
        if self.env.datarepo_version:
            url += "-" + self.env.datarepo_version
        url += ".tar"       # weird, even without this, downloads work

        path = self.generator.bld.bldnode.make_node("tests")
        path = path.abspath()

        debug(f'datarepo: getting {url}')
        web = urlopen(url, context=ssl_context)
        if web.getcode() != 200:
            return web.getcode()

        # tarfile needs seekable file object
        buf = io.BytesIO(web.read())
        buf.seek(0)

        tfname=url.split("/")[-1]
        debug(f'datarepo: unpacking {tfname} in {path}')
        try:
            tf = tarfile.open(fileobj=buf)
            tf.extractall(path)
        except tarfile.TarError as err:
            error(f"datarepo: {err}")
            return -1

        self.outputs[0].write(url)

from waflib.TaskGen import feature
@feature('datarepo')
def create_datarepo_tasks(tgen):
    dlt = tgen.create_task("datarepo_provision")

    dlt.env.datarepo_url = tgen.env.TEST_DATA_URL
    dlt.env.datarepo_category = cat = getattr(tgen, 'category', 'input')
    dlt.env.datarepo_version = ver = getattr(tgen, 'version', '')
    
    if ver:
        urlf = f'{cat}-{ver}.url'
    else:
        urlf = f'{cat}.url'

    tdir = tgen.bld.bldnode.make_node("tests")
    tdir.mkdir()
    dlt.outputs = [tdir.find_or_declare(urlf)]



# class test_data_repo(Task):
#     '''
#     Populate data directories under tests/ in build dir.
#     '''
#     vars = ['TEST_DATA_VERSIONS', 'TEST_DATA_URL']
#     test_data_dir = "tests"
#     def run(self):
#         debug(f"datarepo: downloading from {url}")
#         # Output directory
#         bld = self.generator.bld
#         tdir = bld.bldnode.make_node(self.test_data_dir)
#         tdir.mkdir()
#         path = tdir.abspath()
#         # The history archives.
#         for ver in versions:
#             debug(f"datarepo: populating history/{ver}")
#             rc = download_and_untar(f'{url}/history-{ver}.tar', path)
#             if not rc:
#                 continue
#             return rc
#         # The input archives.
#         # In future, set TEST_DATA_INPUT for versioned input archives.
#         # Unpacking any version must still populate input/.
#         # Include file extension in case future adds compression.
#         # Archive must nonetheless be a tar file.
#         tdi = getattr(self.env, 'TEST_DATA_INPUT', "input.tar")
#         if not isinstance(tdi, str):
#             tdi = tdi[0]
#         debug(f"datarepo: populating input/")
#         return download_and_untar(f'{url}/{tdi}', path)
        
            
def build(bld):
    '''
    Download and unpack archived data repo files.
    '''
    if not bld.env.TESTS:
        return

    bld(features="datarepo", name='datarepo_input')
    for ver in bld.env.TEST_DATA_VERSIONS:
        name = 'datarepo_history_' + ver.replace('.','_')
        bld(features="datarepo", name=name,
            category='history', version=ver)

    # tdir = bld.bldnode.make_node(bld.env.TEST_DATA_DIR)
    # tdir.mkdir()

    # tsk = seed_datarepo_urls(env=bld.env)
    # bld.add_to_group(tsk)
    
    # bld(rule=download_and_untar, name="datarepo_input",
    #     source = tdir.find_or_declare('input.url'))
    # for ver in bld.env.TEST_DATA_VERSIONS:
    #     name = 'datarepo_history_' + ver.replace(".","_")
    #     bld(rule=download_and_untar, name=name,
    #         source = tdir.find_or_declare(f'history-{ver}.url'))



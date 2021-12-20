from pymol import cmd

def parrent(dir):
  return os.path.abspath(os.path.join(dir, os.path.pardir))

scripts = '/home/miller96/Software/dowserplusplus/dpp_original/dowserplusplus/files/result/with_Q/caver_output/1/pymol' + '/'
home = parrent(scripts) + '/'

def exists(name):
	for n in cmd.get_names("all"):
		if n == name:
			return True

def load_safely(file, name):
	if exists(name):
		cmd.delete(name)
	cmd.load(file, name)


view = cmd.get_view()

id = 'cavity_glu_prot_Q10_EM_100kJ_water'
cmd.delete(id + '_*_*')
cmd.delete(id + '_origins')
cmd.delete(id + '_v_origins')

filename = scripts + '/modules/rgb.py'
exec(compile(open(filename, "rb").read(), filename, 'exec'))

color = 1
cluster_dir = home + "data/clusters_timeless"
if os.path.exists(cluster_dir):
	list = os.listdir(cluster_dir)
	list.sort()
	for fn in list:
		name = id + '_' + fn.replace('tun_cl_','t')
		suffix = name[-4:]
		if '.pdb' == suffix or '.ent' == suffix:
			name = name[:-4]
		load_safely(home + 'data/clusters_timeless/' + fn, name)
		cmd.alter(name, 'vdw=b')
		cmd.hide('everything', name)
		cmd.show('spheres', name)
		cmd.color('caver' + str(color), name)
		if color < 1000:
			color += 1

no = id + '_origins'
nvo = id + '_v_origins'
load_safely(home + 'data/origins.pdb', no)
load_safely(home + 'data/v_origins.pdb', nvo)
cmd.show('nb_spheres', no)
cmd.show('nb_spheres', nvo)

cmd.set_view(view)

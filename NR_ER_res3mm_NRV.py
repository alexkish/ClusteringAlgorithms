from ROOT import *
from math import sqrt, pi, fabs
from array import array
from operator import add
import sys

#import psyco
#psyco.full()

#import pdb
#pdb.set_trace()

#set the detector dimensions
sensitive_volume_radius = 154.0 # mm
sensitive_volume_length = 306.5 # mm
cathode_pmt_distance = 14.5 # mm
Zoffset = -21.5

top_veto_low_z = 68.5 # mm
top_veto_radius = 205.0 # mm

bottom_veto_high_z = -378.0 # mm
bottom_veto_radius = 205.0 # mm

#DATA PATH
#schroedinger
#root_files_directory = '/home/physik/alexkish/data/Xe100_AmBe-HP_5e5/'
#idark05x
#root_files_directory = '/Users/alexkish/data/MC/splitChains/'
root_files_directory = '/Users/alexkish/data/MC/neutrons/rock/'

root_file = sys.argv[1]

filename = root_files_directory + root_file
file = TFile(filename, 'read')
tree = file.Get('t1')

nb_entries = tree.GetEntries()

newfilename = root_files_directory
newfile = TFile(root_files_directory + root_file[:-5] + '-NR-ER-res3mm_NRV.root', 'recreate')
newtree = TTree('t2', '')

#che eto takoe??
max_size = 1000

eventid		= array('i', [0])
etot_NR		= array('f', [0.])
etot_ER		= array('f', [0.])
ns_NR		= array('i', [0])
ns_ER		= array('i', [0])
ed_NR 		= array('f', max_size*[0.])
ed_ER 		= array('f', max_size*[0.])
xp_NR 		= array('f', max_size*[0.])
xp_ER 		= array('f', max_size*[0.])
yp_NR 		= array('f', max_size*[0.])
yp_ER 		= array('f', max_size*[0.])
zp_NR 		= array('f', max_size*[0.])
zp_ER 		= array('f', max_size*[0.])
etotv_NR 	= array('f', [0.])
etotv_ER 	= array('f', [0.])
nsv_NR 		= array('i', [0])
nsv_ER 		= array('i', [0])
edv_NR 		= array('f', max_size*[0.])
edv_ER 		= array('f', max_size*[0.])
xpv_NR 		= array('f', max_size*[0.])
xpv_ER 		= array('f', max_size*[0.])
ypv_NR 		= array('f', max_size*[0.])
ypv_ER 		= array('f', max_size*[0.])
zpv_NR 		= array('f', max_size*[0.])
zpv_ER 		= array('f', max_size*[0.])
xp_pri 		= array('f', [0.])
yp_pri 		= array('f', [0.])
zp_pri 		= array('f', [0.])

newtree.Branch('eventid', 	eventid,	'eventid/I')
newtree.Branch('etot_NR', 	etot_NR, 	'etot_NR/F')
newtree.Branch('etot_ER', 	etot_ER, 	'etot_ER/F')
newtree.Branch('ns_NR', 	ns_NR, 		'ns_NR/I')
newtree.Branch('ns_ER', 	ns_ER, 		'ns_ER/I')
newtree.Branch('ed_NR', 	ed_NR, 		'ed_NR[ns_NR]/F')
newtree.Branch('ed_ER', 	ed_ER, 		'ed_ER[ns_ER]/F')
newtree.Branch('xp_NR', 	xp_NR, 		'xp_NR[ns_NR]/F')
newtree.Branch('xp_ER', 	xp_ER, 		'xp_ER[ns_ER]/F')
newtree.Branch('yp_NR', 	yp_NR, 		'yp_NR[ns_NR]/F')
newtree.Branch('yp_ER', 	yp_ER, 		'yp_ER[ns_ER]/F')
newtree.Branch('zp_NR', 	zp_NR, 		'zp_NR[ns_NR]/F')
newtree.Branch('zp_ER', 	zp_ER, 		'zp_ER[ns_ER]/F')
newtree.Branch('etotv_NR', 	etotv_NR, 	'etotv_NR/F')
newtree.Branch('etotv_ER', 	etotv_ER, 	'etotv_ER/F')
newtree.Branch('nsv_NR', 	nsv_NR, 	'nsv_NR/I')
newtree.Branch('nsv_ER', 	nsv_ER, 	'nsv_ER/I')
newtree.Branch('edv_NR', 	edv_NR, 	'edv_NR[nsv_NR]/F')
newtree.Branch('edv_ER', 	edv_ER, 	'edv_ER[nsv_ER]/F')
newtree.Branch('xpv_NR', 	xpv_NR, 	'xpv_NR[nsv_NR]/F')
newtree.Branch('xpv_ER', 	xpv_ER, 	'xpv_ER[nsv_ER]/F')
newtree.Branch('ypv_NR', 	ypv_NR, 	'ypv_NR[nsv_NR]/F')
newtree.Branch('ypv_ER', 	ypv_ER, 	'ypv_ER[nsv_ER]/F')
newtree.Branch('zpv_NR', 	zpv_NR, 	'zpv_NR[nsv_NR]/F')
newtree.Branch('zpv_ER', 	zpv_ER, 	'zpv_ER[nsv_ER]/F')
newtree.Branch('xp_pri', 	xp_pri, 	'xp_pri/F')
newtree.Branch('yp_pri', 	yp_pri, 	'yp_pri/F')
newtree.Branch('zp_pri', 	zp_pri, 	'zp_pri/F')

for entry in range(nb_entries):
#for entry in range(1000):
	tree.GetEntry(entry)

	if entry % 1000 == 0:
		print '(%d/%d)' % (entry, nb_entries)
	
	tree_xp, tree_yp, tree_zp = list(tree.xp), list(tree.yp), list(tree.zp)
	tree_ed = list(tree.ed)
	tree_type = list(tree.type)
	tree_creaproc = list(tree.creaproc)

	# build a list of energy depositions [[e0, [x0, y0, z0]], ...] for all steps, sorted in z
	all_steps = [[tree_ed[step], [tree_xp[step], tree_yp[step], tree_zp[step]], tree_type[step], tree_creaproc[step]] for step in range(tree.nsteps)]
	all_steps.sort(lambda p1, p2: int(p1[1][2]-p2[1][2]))

	# select only steps in the sensitive volume and with energy deposited
	filter_er = lambda p: p[2]=='gamma' or p[2]=='e-' or p[2]=='e+' or p[2]=='mu-' or p[2]=='mu+'
	filter_nr = lambda p: p[2][:2]=='Xe'

	filter_ed = lambda p: p[0] > 0.
	filter_zp = lambda p: p[1][2] < Zoffset and p[1][2] > -sensitive_volume_length+Zoffset
	filter_rp = lambda p: sqrt(p[1][0]**2 + p[1][1]**2) < sensitive_volume_radius

	er_steps = filter(filter_er, all_steps)
	nr_steps = filter(filter_nr, all_steps)

	nz_steps = filter(filter_ed, all_steps)

	nzer_steps = filter(filter_ed, er_steps)
	nznr_steps = filter(filter_ed, nr_steps)

	sver_steps = filter(filter_zp, nzer_steps)
	sver_steps = filter(filter_rp, sver_steps)

	svnr_steps = filter(filter_zp, nznr_steps)
	svnr_steps = filter(filter_rp, svnr_steps)

	# compute energy deposition in the veto
	filter_side_veto = lambda p: sqrt(p[1][0]**2 + p[1][1]**2) > sensitive_volume_radius
	filter_top_veto = lambda p: sqrt(p[1][0]**2 + p[1][1]**2) < top_veto_radius and p[1][2] > top_veto_low_z
	filter_bottom_veto = lambda p: sqrt(p[1][0]**2 + p[1][1]**2) < bottom_veto_radius and p[1][2] < bottom_veto_high_z

	side_vetonr_steps = filter(filter_side_veto, nznr_steps)
	top_vetonr_steps = filter(filter_top_veto, nznr_steps)
	bottom_vetonr_steps = filter(filter_bottom_veto, nznr_steps)

	side_vetoer_steps = filter(filter_side_veto, nzer_steps)
	top_vetoer_steps = filter(filter_top_veto, nzer_steps)
	bottom_vetoer_steps = filter(filter_bottom_veto, nzer_steps)

	vetonr_steps = side_vetonr_steps
	vetonr_steps += [step for step in top_vetonr_steps if step not in vetonr_steps]
	vetonr_steps += [step for step in bottom_vetonr_steps if step not in vetonr_steps]

	vetoer_steps = side_vetoer_steps
	vetoer_steps += [step for step in top_vetoer_steps if step not in vetoer_steps]
	vetoer_steps += [step for step in bottom_vetoer_steps if step not in vetoer_steps]

	############################################################
	# compute the number of scatters in the sensitive volume, NR
	if len(svnr_steps) > 0:
		# start with 1 scatter
		nbnr_s = 1

		snr_etot = svnr_steps[0][0]
			
		snr_ed = svnr_steps[0][0]
		snr_xp = svnr_steps[0][1][0]*svnr_steps[0][0]
		snr_yp = svnr_steps[0][1][1]*svnr_steps[0][0]
		snr_zp = svnr_steps[0][1][2]*svnr_steps[0][0]

		# for all steps starting with the second one
		for step in range(1, len(svnr_steps)):
			stepnr_ed = svnr_steps[step][0]
			stepnr_xp, stepnr_yp, stepnr_zp = svnr_steps[step][1][0], svnr_steps[step][1][1], svnr_steps[step][1][2]

			# if the step if separated enough
			if fabs(stepnr_zp - svnr_steps[step-1][1][2]) > 3.0:
				# then this is a new scatter

				# save the energy deposited and the weighted position of the previous scatter
				ed_NR[nbnr_s-1] = snr_ed
				xp_NR[nbnr_s-1] = snr_xp/snr_ed
				yp_NR[nbnr_s-1] = snr_yp/snr_ed
				zp_NR[nbnr_s-1] = snr_zp/snr_ed

				snr_etot += stepnr_ed

				# start a new scatter
				nbnr_s += 1
				snr_ed = stepnr_ed
				snr_xp, snr_yp, snr_zp = stepnr_xp*stepnr_ed, stepnr_yp*stepnr_ed, stepnr_zp*stepnr_ed

			else:
				# this is not a new scatter, keep adding in the current
				snr_etot += stepnr_ed
				snr_ed += stepnr_ed
				snr_xp += stepnr_xp*stepnr_ed
				snr_yp += stepnr_yp*stepnr_ed
				snr_zp += stepnr_zp*stepnr_ed

		# save the last scatter (or the first if we have only one!)
		ed_NR[nbnr_s-1] = snr_ed
		xp_NR[nbnr_s-1] = snr_xp/snr_ed
		yp_NR[nbnr_s-1] = snr_yp/snr_ed
		zp_NR[nbnr_s-1] = snr_zp/snr_ed
	else:
		nbnr_s = 0


	############################################################
	# compute the number of scatters in the sensitive volume, ER
	if len(sver_steps) > 0:
		# start with 1 scatter
		nber_s = 1

		ser_etot = sver_steps[0][0]
			
		ser_ed = sver_steps[0][0]
		ser_xp = sver_steps[0][1][0]*sver_steps[0][0]
		ser_yp = sver_steps[0][1][1]*sver_steps[0][0]
		ser_zp = sver_steps[0][1][2]*sver_steps[0][0]

		# for all steps starting with the second one
		for step in range(1, len(sver_steps)):
			steper_ed = sver_steps[step][0]
			steper_xp, steper_yp, steper_zp = sver_steps[step][1][0], sver_steps[step][1][1], sver_steps[step][1][2]

			# if the step if separated enough
			if fabs(steper_zp - sver_steps[step-1][1][2]) > 3.0:
				# then this is a new scatter

				# save the energy deposited and the weighted position of the previous scatter
				ed_ER[nber_s-1] = ser_ed
				xp_ER[nber_s-1] = ser_xp/ser_ed
				yp_ER[nber_s-1] = ser_yp/ser_ed
				zp_ER[nber_s-1] = ser_zp/ser_ed

				ser_etot += steper_ed

				# start a new scatter
				nber_s += 1
				ser_ed = steper_ed
				ser_xp, ser_yp, ser_zp = steper_xp*steper_ed, steper_yp*steper_ed, steper_zp*steper_ed

			else:
				# this is not a new scatter, keep adding in the current
				ser_etot += steper_ed
				ser_ed += steper_ed
				ser_xp += steper_xp*steper_ed
				ser_yp += steper_yp*steper_ed
				ser_zp += steper_zp*steper_ed

		# save the last scatter (or the first if we have only one!)
		ed_ER[nber_s-1] = ser_ed
		xp_ER[nber_s-1] = ser_xp/ser_ed
		yp_ER[nber_s-1] = ser_yp/ser_ed
		zp_ER[nber_s-1] = ser_zp/ser_ed
	else:
		nber_s = 0

	################################################
	# compute the number of scatters in the veto, NR
	if len(vetonr_steps) > 0:
		# start with 1 scatter
		nbnr_vs = 1

		vsnr_etot = vetonr_steps[0][0]
			
		vsnr_ed = vetonr_steps[0][0]
		vsnr_xp = vetonr_steps[0][1][0]*vetonr_steps[0][0]
		vsnr_yp = vetonr_steps[0][1][1]*vetonr_steps[0][0]
		vsnr_zp = vetonr_steps[0][1][2]*vetonr_steps[0][0]

		# for all steps starting with the second one
		for step in range(1, len(vetonr_steps)):
			stepnr_ed = vetonr_steps[step][0]
			stepnr_xp, stepnr_yp, stepnr_zp = vetonr_steps[step][1][0], vetonr_steps[step][1][1], vetonr_steps[step][1][2]

			# if the step if separated enough
			if fabs(sqrt((stepnr_xp - vetonr_steps[step-1][1][0])**2 + (stepnr_yp - vetonr_steps[step-1][1][1])**2 + (stepnr_zp - vetonr_steps[step-1][1][2])**2)) > 30.0:
				# then this is a new scatter

				# save the energy deposited and the weighted position of the previous scatter
				edv_NR[nbnr_vs-1] = vsnr_ed
				xpv_NR[nbnr_vs-1] = vsnr_xp/vsnr_ed
				ypv_NR[nbnr_vs-1] = vsnr_yp/vsnr_ed
				zpv_NR[nbnr_vs-1] = vsnr_zp/vsnr_ed

				vsnr_etot += stepnr_ed

				# start a new scatter
				nbnr_vs += 1
				vsnr_ed = stepnr_ed
				vsnr_xp, vsnr_yp, vsnr_zp = stepnr_xp*stepnr_ed, stepnr_yp*stepnr_ed, stepnr_zp*stepnr_ed

			else:
				# this is not a new scatter, keep adding in the current
				vsnr_etot += stepnr_ed
				vsnr_ed += stepnr_ed
				vsnr_xp += stepnr_xp*stepnr_ed
				vsnr_yp += stepnr_yp*stepnr_ed
				vsnr_zp += stepnr_zp*stepnr_ed

		# save the last scatter (or the first if we have only one!)
		edv_NR[nbnr_vs-1] = vsnr_ed
		xpv_NR[nbnr_vs-1] = vsnr_xp/vsnr_ed
		ypv_NR[nbnr_vs-1] = vsnr_yp/vsnr_ed
		zpv_NR[nbnr_vs-1] = vsnr_zp/vsnr_ed
	else:
		nbnr_vs = 0
		vsnr_etot = 0

	################################################
	# compute the number of scatters in the veto, ER
	if len(vetoer_steps) > 0:
		# start with 1 scatter
		nber_vs = 1

		vser_etot = vetoer_steps[0][0]
			
		vser_ed = vetoer_steps[0][0]
		vser_xp = vetoer_steps[0][1][0]*vetoer_steps[0][0]
		vser_yp = vetoer_steps[0][1][1]*vetoer_steps[0][0]
		vser_zp = vetoer_steps[0][1][2]*vetoer_steps[0][0]

		# for all steps starting with the second one
		for step in range(1, len(vetoer_steps)):
			steper_ed = vetoer_steps[step][0]
			steper_xp, steper_yp, steper_zp = vetoer_steps[step][1][0], vetoer_steps[step][1][1], vetoer_steps[step][1][2]

			# if the step if separated enough
			if fabs(sqrt((steper_xp - vetoer_steps[step-1][1][0])**2 + (steper_yp - vetoer_steps[step-1][1][1])**2 + (steper_zp - vetoer_steps[step-1][1][2])**2)) > 30.0:
				# then this is a new scatter

				# save the energy deposited and the weighted position of the previous scatter
				edv_ER[nber_vs-1] = vser_ed
				xpv_ER[nber_vs-1] = vser_xp/vser_ed
				ypv_ER[nber_vs-1] = vser_yp/vser_ed
				zpv_ER[nber_vs-1] = vser_zp/vser_ed

				vser_etot += steper_ed

				# start a new scatter
				nber_vs += 1
				vser_ed = steper_ed
				vser_xp, vser_yp, vser_zp = steper_xp*steper_ed, steper_yp*steper_ed, steper_zp*steper_ed

			else:
				# this is not a new scatter, keep adding in the current
				vser_etot += steper_ed
				vser_ed += steper_ed
				vser_xp += steper_xp*steper_ed
				vser_yp += steper_yp*steper_ed
				vser_zp += steper_zp*steper_ed

		# save the last scatter (or the first if we have only one!)
		edv_ER[nber_vs-1] = vser_ed
		xpv_ER[nber_vs-1] = vser_xp/vser_ed
		ypv_ER[nber_vs-1] = vser_yp/vser_ed
		zpv_ER[nber_vs-1] = vser_zp/vser_ed
	else:
		nber_vs = 0
		vser_etot = 0


	########################
	# save event information
	if nbnr_s>0 and nber_s>0:
		eventid[0] = tree.eventid
		etot_NR[0] = snr_etot
		etot_ER[0] = ser_etot
		ns_NR[0] = nbnr_s
		ns_ER[0] = nber_s
		# veto
		etotv_NR[0] = vsnr_etot 
		etotv_ER[0] = vser_etot 
		nsv_NR[0] = nbnr_vs
		nsv_ER[0] = nber_vs
		# primaries
		xp_pri[0] = tree.xp_pri
		yp_pri[0] = tree.yp_pri
		zp_pri[0] = tree.zp_pri

		# fill the tree
		newtree.Fill()

	if nbnr_s==0 and nber_s>0:
		eventid[0] = tree.eventid
		etot_NR[0] = 0
		etot_ER[0] = ser_etot
		ns_NR[0] = 0
		ns_ER[0] = nber_s

		etotv_NR[0] = vsnr_etot
		etotv_ER[0] = vser_etot
		nsv_NR[0] = nbnr_vs
		nsv_ER[0] = nber_vs
		newtree.Fill()


	if nbnr_s>0 and nber_s==0:
		# save event information
		eventid[0] = tree.eventid
		etot_NR[0] = snr_etot
		etot_ER[0] = 0
		ns_NR[0] = nbnr_s
		ns_ER[0] = 0

		etotv_NR[0] = vsnr_etot 
		etotv_ER[0] = vser_etot 
		nsv_NR[0] = nbnr_vs
		nsv_ER[0] = nber_vs

		xp_pri[0] = tree.xp_pri
		yp_pri[0] = tree.yp_pri
		zp_pri[0] = tree.zp_pri

		# fill the tree
		newtree.Fill()

file.Close()

newfile.Write()
newfile.Close()


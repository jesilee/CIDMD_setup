############################################################################
# Jesi Lee                                                                 #
# Jan 2022		                                                   #
#                                                                          #
# This is the master script that creates the subfolderse cid/calcs         #
#									   #
# Please set/change the variable options				   #
#					                                   #
# dependencies: jesi_SetupRanCID3.py					   #
#		jesi_PrepRanCID.py   					   #
#		xyz file                                                   #
#                                                                          #
#                                                                          #
############################################################################


my_opt_xyz=BF_case1_optim.xyz

python3 CIDMD_setup.py --n_points 200 --v_scale 8 --cutoff 0.5 $my_opt_xyz

python3 CIDMD_preprun.py $my_opt_xyz

mkdir cid 

mv calcs/ cid/

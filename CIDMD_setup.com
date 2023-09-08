############################################################################
# Jesi Lee                                                                 
# Jan 2022		                                                   
#   
# Please set/change the variable options				   
#					                                   
#                                                                          
############################################################################


my_opt_xyz=PS_case1_optim.xyz

python3 CIDMD_setup.py --n_points 200 --v_scale 8 --cutoff 0.5 $my_opt_xyz

python3 CIDMD_preprun.py $my_opt_xyz

mkdir cid 

mv calcs/ cid/

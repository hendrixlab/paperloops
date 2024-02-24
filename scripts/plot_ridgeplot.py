import matplotlib.pyplot as plt
import pandas as pd
import joypy

pl_pts_fwd = pd.read_csv('pts_fwd_file.txt',sep='\t')
pl_pts_fwd = pl_pts_fwd.loc[:, pl_pts_fwd.columns != 'idx']

tf3_pts_fwd = pd.read_csv('../toyfold3/runs/sampled_pts_fwd.txt',sep='\t')
tf3_pts_fwd = tf3_pts_fwd.loc[:, tf3_pts_fwd.columns != 'idx']

df = pd.concat([pl_pts_fwd, tf3_pts_fwd], axis=1, join='inner')

df = pd.DataFrame(df, columns=['pl_x','tf3_x','pl_y','tf3_y','pl_z','tf3_z','pl_vx','tf3_vx','pl_vy','tf3_vy','pl_vz','tf3_vz'])


plt.figure()

fig2, ax2 = joypy.joyplot(df)

plt.savefig('pts_fwd_compared.pdf')
plt.clf()

pl_pts_rev = pd.read_csv('pts_rev_file.txt',sep='\t')
pl_pts_rev = pl_pts_rev.loc[:, pl_pts_rev.columns != 'idx']

tf3_pts_rev = pd.read_csv('../toyfold3/runs/sampled_pts_rev.txt',sep='\t')
tf3_pts_rev = tf3_pts_rev.loc[:, tf3_pts_rev.columns != 'idx']

df = pd.concat([pl_pts_rev, tf3_pts_rev], axis=1, join='inner')

df = pd.DataFrame(df, columns=['pl_x','tf3_x','pl_y','tf3_y','pl_z','tf3_z','pl_vx','tf3_vx','pl_vy','tf3_vy','pl_vz','tf3_vz'])

plt.figure()

fig2, ax2 = joypy.joyplot(df)

plt.savefig('pts_rev_compared.pdf')

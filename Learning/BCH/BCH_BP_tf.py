import tensorflow as tf
import numpy as np
import scipy.sparse as sps

n=15
m=4 #Number of check nodes   

H=np.array([[1,0,0,0,1,1,1,0,0,0,0,1,1,1,1],
	        [0,1,0,0,1,0,0,1,1,0,1,0,1,1,1],
		    [0,0,1,0,0,1,0,1,0,1,1,1,0,1,1],
			[0,0,0,1,0,0,1,0,1,1,1,1,1,0,1]])

(m1,n1)=H.shape
assert n1==n and m1==m
del m1, n1



x=t.placeholder(tf.float32,shape=[None,n])
y_=tf.placeholder(tf.boolean,shape=[None,n])




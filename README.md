# Condition-Driven_Process_Monitoring
Source code of Condition-Driven Data Analytics and Monitoring for Wide-Range Nonstationary and Transient Continuous Processes.   
The details of model can be found in    
 [C. Zhao, J. Chen and H. Jing, "Condition-Driven Data Analytics and Monitoring for Wide-Range Nonstationary and Transient Continuous Processes," in IEEE Transactions on Automation Science and Engineering, vol. 18, no. 4, pp. 1563-1574, Oct. 2021.](https://ieeexplore.ieee.org/abstract/document/9158352)

#### Example:  
Please see 'demo.m' for how to use this model.
The data used in this paper is not allowed to be shared. You should prepare the data and change the parameters accordingly before run 'demo.m'.

#### Note:
* We use monitoring methods SFA/PCA as base models and monitoring statistics as merge indicators to capture process characteristics and divide data into different modes. You can change the base model and statistics according to your needs. If so, you should prepare your own class based on the 'base_model/SFA_class.m'.
* If you want to segment and rearrange data for other tasks rather than monitoring tasks, you simply need to adjust the methods and indicators accordingly. For example, if segmented for regression tasks, SFA and PCA can be replaced with regression methods such as LR, and the merging indicator can be set to regression errors.
* The 'demo.m' is an example showing how to use divided data for monitoring. You can use divided data for other purposes. If so, you should replace the 'utils/monitoring.m' function with your own function. 

#### All rights reserved, citing the following paper is required for reference:   
[1] C. Zhao, J. Chen and H. Jing, "Condition-Driven Data Analytics and Monitoring for Wide-Range Nonstationary and Transient Continuous Processes," in IEEE Transactions on Automation Science and Engineering, vol. 18, no. 4, pp. 1563-1574, Oct. 2021.  

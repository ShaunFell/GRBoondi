## GRBoondi 2024
## Copyright 2024, Shaun Fell
## Please refer to LICENSE in GRBoondi's root directory


## The engine launcher sets up the parallel engine to generate plots 
## on high performance computing clusters that utilize job scheduling (e.g. slurm).
## If parallelization is not desired, users simply have to set 'use_parallel=false' 
## in their parameter file. The engine will instead be launched on the local machine 
## that the script is executed on. 


def setup_engine(config):
	"""Launch the engine for generating the plots. Uses job scheduling on high performance computing clusters. Uses the local machine if 'useparallel' is false.

	Args:
		config (configparser.ConfigParser): instance of a ConfigParser class that holds the users parameters
	"""
	useparallel = config["EngineConfig"].getboolean("use_parallel", 0)
	
	#if useparallel is disabled, default engine is launched
	if useparallel:
		host = config["EngineConfig"].get("host", "localhost")
		print("Running in parallel on host: ", host)

		#Parse config options for engine
		num_procs = config["EngineConfig"].get("number_processes", 1)
		num_nodes = config["EngineConfig"].get("number_nodes", 1)
		partition_name = config["EngineConfig"].get("partition")
		job_cmd = config["EngineConfig"].get("job_cmd", "srun")
		time_limit = config["EngineConfig"].get("time_limit", "01:00:00")
		job_name = config["EngineConfig"].get("job_name", "3DPlot")
		add_sub_args = config["EngineConfig"].get("additional_launch_args", "")
		use_gpus = config["EngineConfig"].getboolean("use_gpus", 0)
		ngpus_per_node = config["EngineConfig"].get("ngpus_per_node", "1")

		#create argument tuple
		arg = ("-l", job_cmd, "-n", job_name, "-p", partition_name, "-np", num_procs, "-nn", num_nodes, "-t", time_limit, "-la", add_sub_args)

		#if gpus requested, add appropriate arguments to list
		if use_gpus:
			slurm_gpu_submission = "--gres=gpu:{0}".format(ngpus_per_node)
			add_sub_args += " {0}".format(slurm_gpu_submission)
			arg = arg[:-2] + ("-la", add_sub_args,)
			arg = arg + ("-hw_accel",)
			arg = arg + ("-n-gpus-per-node", ngpus_per_node,)
		

		#Submit the job
		openengine_status = OpenComputeEngine(host, arg)

        #if job submission fails and user desires a parallelized execution, then we must exit
		if not openengine_status:
			print("Job submission failed. Exiting...")
			sys.exit(1)


import urllib
import string
import os
import sys
import LaunchOnCondor

FarmDirectory = "condor"
JobName = "cuba"
LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)

path = "/home/fynu/swertz/tests_MEM/MEMcpp/"

start_evt = 0
end_evt = -1
evt_per_job = 50
max_evt= 2000
i = 0

while start_evt < max_evt:
	end_evt += evt_per_job
	LaunchOnCondor.SendCluster_Push(["BASH", path + "/condor.sh", path + "/data/ttbar_weighted_noTF_isr0_pdfMtop_noMCoPerms_2000evt.root", "ttbar_noTF_isr0_pdfMtop_noMCoPerms_as013_cuba_Mersenne_noSmooth_50000_" + str(i) + ".root", start_evt, end_evt])
	#LaunchOnCondor.SendCluster_Push(["BASH", path + "/condor.sh", path + "/results/ttbar_noTF_isr0_pdfMtop_noMCoPerms_as013_cuba_Mersenne_noSmooth_50000.root", "ttbar_noTF_isr0_pdfMtop_noMCoPerms_as013_cuba_Mersenne_noSmooth_50000_bis_" + str(i) + ".root", start_evt, end_evt])
	start_evt = end_evt+1
	i += 1

LaunchOnCondor.SendCluster_Submit()

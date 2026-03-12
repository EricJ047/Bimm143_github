## Basic Unix

Some important file system commands include 

pwd: print working directory 
ls: list files and folders
cd: change directory 
mdkir: make a new directory
rm: delete files and directories*
nano: a very basic text editor that is always available
less: to view/read text files page by page(pager program)

My AWS instance:
ssh -i ~/downloads/BIMM143_YJ.pem ubuntu@ec2-35-95-26-97.us-west-2.compute.amazonaws.com

To copy from my AWS instance:
scp -i ~/downloads/BIMM143_YJ.pem ubuntu@ec2-35-95-26-97.us-west-2.compute.amazonaws.com:~/work/results.tsv .

Connect:
ssh -i "BIMM143_YJ.pem" ubuntu@ec2-35-85-57-248.us-west-2.compute.amazonaws.com
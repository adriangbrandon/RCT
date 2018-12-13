import json

file = open("rct_index_40_1024.json")
data = json.load(file)

#Log Object
print data["children"][9]["name"]
log_object_data = data["children"][9]["children"][0]
#--> Disappearances

for c in log_object_data["children"]:
	print c["name"] + "\t" + c["size"]
print "\n"


#Log Reference
print data["children"][10]["name"]
log_reference_data = data["children"][10]
#--> 
for c in log_reference_data["children"]:
	print c["name"] + "\t" + c["size"]
print "\n"


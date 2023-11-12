def readpdbfile(input):
        residue = []
        #molnumber = []
        #atomname = []
        atomtype = []
        groupname = []

        with open(input) as pdbfile:
                for line in pdbfile:
                        if line[:4] == 'ATOM' or line[:6] == "HETATM":
                                #print(line) #DEBUG
                                # Split the line
                                splitted_line = [line[:6], line[6:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38], line[38:46], line[46:54]]
                                #print(splitted_line) #DEBUG
                                # To format again the pdb file with the fields extracted
                                #print("%-6s%5s %4s %3s %s%4s    %8s%8s%8s\n"%tuple(splitted_line)) #DEBUG
                                residue.append(splitted_line[3].strip())
                                atomtype.append(splitted_line[2].strip())
                                #atomname.append(splitted_line[11].strip())
                                #molnumber.append(int(splitted_line[5])) 
                                groupname.append(splitted_line[4].strip())
                                #print(splitted_line[5]) #DEBUG
                uniqueres = []
                for res in residue:
                    if res not in uniqueres:
                        uniqueres.append(res)
                resnum = len(uniqueres)
        return residue, atomtype, groupname, resnum, uniqueres


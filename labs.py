def loadBlosum50():
    with open('blosum50.txt') as f:
        content = f.readlines()
    content = [x.strip() for x in content]

    blosumCosts = []
    for x in content:
        lineArray = x.split(' ')
        lineArray = [int(x) for x in lineArray if x]
        blosumCosts.append(lineArray)
    return blosumCosts

def subCost(a, b):
    return 5

def nwAlgo(a, b):
    insertionCost = -8
    costMatrix = []
        
    for i in range(len(a)):
        costMatrix.append([])
        for j in range(len(b)):
            vals = []
            if(i > 0 and j > 0):
                vals.append(costMatrix[i-1][j-1] + subCost(a[i],b[j]))

            if(i > 0):
                vals.append(costMatrix[i-1][j] + insertionCost)

            if(j > 0):
                vals.append(costMatrix[i][j-1] + insertionCost)

            if not vals:
                vals.append(0)
            
            costMatrix[i].append(max(vals))

    return costMatrix

matrix = nwAlgo("abcde", "abcde")
for row in matrix:
    print(row)

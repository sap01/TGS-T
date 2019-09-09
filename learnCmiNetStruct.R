#------------------
#CMI Implemetation:
#------------------ 
learnCmiNetStruct <- function(di.net.adj.matrix, input.data, num.nodes, beta)
{
  for(i in 1:num.nodes) # for each node i
  {
    adj_list_1 <- c()
    adj_list_2 <- c()
    
    for(j in 1:num.nodes) # for each node j
    {
      # if i is a parent of j in 'final_net'
      if(di.net.adj.matrix[i,j] == 1)
      {
        for(k in 1:num.nodes) # for each node k
        {
          # if 'k' is a parent of 'i' and a grandparent of 'j'
          # i.e. there exists a directed chain 
          # from 'k' to 'i' to 'j' in 'final_net'
          # ToDo: Also check '(k != i)' because 'k' may be 
          # a parent as well as a grandparent of 'i'
          if(di.net.adj.matrix[k,i] == 1 && k != j)
          {
            # Append 'k' to 'adj_list_1'
            # 'adj_list_1' contains all grandparents of 'j' except 'j'
            adj_list_1 <- c(adj_list_1, k)
          }
        }
        #ToDo: Can not we merge this for loop with the prev for loop?
        for(k in 1:num.nodes) # for each node k
        {
          # if 'k' (!= 'i') is a parent of 'j'
          if(di.net.adj.matrix[k,j] == 1 && k != i)
          {
            # 'adj_list_2' contains all parents of 'j' except 'i'
            adj_list_2 <- c(adj_list_2, k)
          }
        }
        
        # 'commom_node' contains all nodes except 'i' that are grandparents
        # as well as parents of 'j' 
        common_node <- intersect(adj_list_1,adj_list_2)
        
        #print(length(common_node)) #---test
        
        if(length(common_node) > 0)
        {
          max_I <- 0
          
          # Store max(CMI('i', 'j' | 'node') for all 'node' in 'commom_node')
          # in max_I
          for (node in common_node)
          {
            # Calc CMI('i', 'j' | 'node')
            # Another potential point of difference from the original LBN implementaiton.
            # Here, "emp" estimator is used to calc CMI. The orig implementaion migh
            # have used a diff estimator.
            I <- computeCmi(input.data[,i],input.data[,j],input.data[,node])
            
            if(I > max_I)
            {
              max_I <- I
            }	
          }
          #2nd order CMI		
          if(length(common_node) > 1)
          {
            # All combination of the nodes in 'commom_node' taken 2 at a time.
            # In 'comb2', row = memeber index, col = combination index.
            comb_2 <- combn(common_node, 2)
            
            # ToDo: 'colunm in 1:ncol(comb_2)'
            # ToDo: Correct spelling 'colunm'
            for (column in 1:ncol(comb_2))
            {
              # ToDo: I <- condinformation(dat[,i],dat[,j],dat[, comb_2[, colunm]],method="emp")
              #I <- condinformation(dat[,i],dat[,j],dat[,column],method="emp")
              I <- computeCmi(input.data[,i],input.data[,j],input.data[, comb_2[, column]])
              #	print(I)  #---test
              if(I > max_I)
              {
                max_I <- I
              }	
            }						
            
          }
          
          # ToDo: Differentiate alpha and beta
          # ToDo: Instead of recording 'max_I', use 'break' construct to break out
          # of a loop when 'I < alpha' 
          if(max_I < beta)
          {
            # Remove the directed edge from 'i' to 'j'
            di.net.adj.matrix[i,j] <- 0
          }	
        }
      }
    }
  }
  return (di.net.adj.matrix)
}
#########################
#  GF(4) Arithmetic
#########################

	function matrixGFMult(matrixA,matrixB)
	

		n1=length(matrixA[:,1])
		n2=length(matrixB[:,1])

		m1=length(matrixA[1,:])
		m2=length(matrixB[1,:])
		
		matrixRes=Int.(zeros(n1,m2))

		if(m1!=n2)
			println("error during matrix GF multiplication: number of columns ( matrix A ",m1,") different from number of rows (matrix B ",n2,")")
			return
		else

			for i=1:n1
				for j=1:m2
					for k=1:m1
						tmp=multGF(matrixA[i,k],matrixB[k,j])
						matrixRes[i,j]=addGF(matrixRes[i,j],tmp)
					end				
				end
			end


		end

		return matrixRes

	end

	function divGF(symb1,symb2)
		res=0
		if(symb1==0)
			if(symb2==0)
				res=0
			elseif(symb2==1)
				res=0
			elseif(symb2==2)
				res=0
			elseif(symb2==3)
				res=0
			end
		elseif(symb1==1)
			if(symb2==0)
				res=0
			elseif(symb2==1)
				res=1
			elseif(symb2==2)
				res=3
			elseif(symb2==3)
				res=2
			end
		elseif(symb1==2)
			if(symb2==0)
				res=0
			elseif(symb2==1)
				res=2
			elseif(symb2==2)
				res=1
			elseif(symb2==3)
				res=3
			end
		elseif(symb1==3)
			if(symb2==0)
				res=0
			elseif(symb2==1)
				res=3
			elseif(symb2==2)
				res=2
			elseif(symb2==3)
				res=1
			end
		end
		return res
	end


	function addGF(symb1,symb2)
		res=0
		if(symb1==0)
			if(symb2==0)
				res=0
			elseif(symb2==1)
				res=1
			elseif(symb2==2)
				res=2
			elseif(symb2==3)
				res=3
			end
		elseif(symb1==1)
			if(symb2==0)
				res=1
			elseif(symb2==1)
				res=0
			elseif(symb2==2)
				res=3
			elseif(symb2==3)
				res=2
			end
		elseif(symb1==2)
			if(symb2==0)
				res=2
			elseif(symb2==1)
				res=3
			elseif(symb2==2)
				res=0
			elseif(symb2==3)
				res=1
			end
		elseif(symb1==3)
			if(symb2==0)
				res=3
			elseif(symb2==1)
				res=2
			elseif(symb2==2)
				res=1
			elseif(symb2==3)
				res=0
			end
		end
		return res
	end

	function multGF(symb1,symb2)
		res=0
		if(symb1==0)
			if(symb2==0)
				res=0
			elseif(symb2==1)
				res=0
			elseif(symb2==2)
				res=0
			elseif(symb2==3)
				res=0
			end
		elseif(symb1==1)
			if(symb2==0)
				res=0
			elseif(symb2==1)
				res=1
			elseif(symb2==2)
				res=2
			elseif(symb2==3)
				res=3
			end
		elseif(symb1==2)
			if(symb2==0)
				res=0
			elseif(symb2==1)
				res=2
			elseif(symb2==2)
				res=3
			elseif(symb2==3)
				res=1
			end
		elseif(symb1==3)
			if(symb2==0)
				res=0
			elseif(symb2==1)
				res=3
			elseif(symb2==2)
				res=1
			elseif(symb2==3)
				res=2
			end
		end
		return res

	end

#### END GF CALCULATION #####


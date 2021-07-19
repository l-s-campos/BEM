## Início da análise
# include("potencial_interpolado.jl")
using TimerOutputs
to = TimerOutput()
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
nelem = 100  #Numero de elementos
NPX = 5 #pontos internos na direção x
NPY = 5 #pontos internos na direção y
npg = 16    #apenas números pares
## Formatação dos dados ________________________________________________
println("1. Formatando os dados");
@timeit to "dados" dad = format_dad_iga(potencial1diso(nelem),NPX,NPY) # dados

println("2. Montando a matriz A e o vetor b")
@timeit to "padrão" H,G = calc_HeG(dad,npg)  #importante
# Ht,Gt = calc_HeG(dad,1:nc(dad)+ni(dad),1:4*nelem,npg)  #importante
@timeit to "arvore" Tree1,Tree2,block=cluster(dad,η = 1.1,max_elem=30)
# Haca,Baca=Hinterp(dad,Tree1,Tree2,block)
# @timeit to "interp" Ai,bi=Ainterp(dad,Tree1,Tree2,block,ninterp=5)
@timeit to "kmeans" A1,b1=Akmeans(dad,Tree1,Tree2,block,nnucleos=10);
# @timeit to "kmeans2" A2,b2=Akmeans2(dad,Tree1,Tree2,block,nnucleos=9)
show(to)
# Aicheia=BEM.montacheia(Ai,block,Tree1,Tree2,dad,Ht)
A,b = aplicaCDC(H,G,dad) # Calcula a matriz A e o vetor b
println("3. Resolvendo o sistema linear")
# x2,f = gmres(A,b,5,tol=1e-5,maxIter=100,out=0) #GMRES nas matrizes hierarquicas
x = A\b
println("4. Separando fluxo e temperatura")
T,q = separa(dad,x) #importante
# A1,b1 = calc_Aeb(dad,npg)  #importante
# tol =1e-6
# @time x1=A1\b1;
# @time x2,r = gmres(A,b ,3,tol=tol,maxIter=100,out=0);
# @time x3,r =fgmres(A,b ,3,tol=tol,maxIter=100,out=0);
# @time x4,r =bicgstb(A,b ,tol=tol,maxIter=100,out=0);
# @time x5,r =minres(A,b ,rtol=tol,maxIter=100,out=0);
# Ac = BEM.montaAcheia(Ai,block,Tree1,Tree2,dad)#matriz cheia
# Ac1= BEM.montaAcheia(A1,block,Tree1,Tree2,dad)#matriz cheia
# Ac2= BEM.montaAcheia(A2,block,Tree1,Tree2,dad)#matriz cheia
# @show norm(Ac[1:nc(dad),1:nc(dad)]-A[1:nc(dad),1:nc(dad)])
# @show norm(Ac1[1:nc(dad),1:nc(dad)]-A[1:nc(dad),1:nc(dad)])
# @show norm(Ac2[1:nc(dad),1:nc(dad)]-A[1:nc(dad),1:nc(dad)])
# [bi[1:nc(dad)] b1[1:nc(dad)] b2[1:nc(dad)]].-b
# Pre= BEM.montaAcheia(Ai,block,Tree1,Tree2,dad,full=false)#matriz pre dos blocos não admissiveis

# xi,f = gmres(vet->matvec(Ai,vet,block,Tree1,Tree2,dad),bi,5,tol=1e-5,maxIter=100,out=0) #GMRES nas matrizes hierarquicas
x1,f = gmres(vet->matvec(A1,vet,block,Tree1,Tree2,dad),b1,5,tol=1e-5,maxIter=100,out=0) #GMRES nas matrizes hierarquicas
T1,q1 = separa(dad,x1) #importante
T1i=x1[nc(dad)+1:end]
println("5. Calculando nos pontos internos")
Ti = calc_Ti(dad,T,q,npg);

println("6. Gerando mapa de cor")
erro=norm([T-dad.NOS[:,1];Ti-dad.pontos_internos[:,1]])
erroi=norm([T1-dad.NOS[:,1];T1i-dad.pontos_internos[:,1]])

# geo=mostra_geometria(dad)
# mapa=mostra_resultado(dad,[T;Ti])

# @btime matvec(Aaca,baca,block,Tree1,Tree2,dad)
# @btime Ac*baca
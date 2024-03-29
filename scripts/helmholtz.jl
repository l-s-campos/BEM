## Início da análise
using DrWatson
@quickactivate "BEM"
include(scriptsdir("includes.jl"))
nelem = 1  #Numero de elementos
NPX = 2 #pontos internos na direção x
NPY = 2 #pontos internos na direção y
npg = 10    #apenas números pares
## Formatação dos dados ________________________________________________
println("1. Formatando os dados");
dad = format_dad(helm1d(nelem,3),NPX,NPY) # dados
# dad0 = format_dad(potencial1d(nelem),NPX,NPY,0) # dados
# dad = format_dad(placacomfuro(nelem),NPX,NPY) # dados

println("2. Montando a matriz A e o vetor b")
H,G = calc_HeG(dad,10)  #importante
A,b = aplicaCDC(H,G,dad) # Calcula a matriz A e o vetor b
println("3. Resolvendo o sistema linear")
x = A\b
println("4. Separando fluxo e temperatura")
T,q = separa(dad,x) #importante

T1=real(T)
q1=real(q)
	
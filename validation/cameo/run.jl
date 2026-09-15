using MizuRoute

work = length(ARGS)>=1 ? ARGS[1] : "validation/cameo/work"
meta_txt = read(joinpath(work,"meta.json"), String)
function num(key)
    m=match(Regex("\\\""*key*"\\\"\\s*:\\s*([-+0-9.eE]+)"),meta_txt)
    m===nothing && error("missing meta key $key")
    parse(Float64,m.captures[1])
end
nreach=Int(num("nreach")); ntime=Int(num("ntime")); dt=num("dt")
lines=readlines(joinpath(work,"network.csv"))[2:end]
rid=Vector{Int}(undef,nreach); down=similar(rid)
len=zeros(nreach); slope=zeros(nreach); width=zeros(nreach); active=falses(nreach)
for (i,line) in enumerate(lines)
    x=split(line,',')
    rid[i]=parse(Int,x[1]); down[i]=parse(Int,x[2]); len[i]=parse(Float64,x[3])
    slope[i]=parse(Float64,x[4]); width[i]=parse(Float64,x[5]); active[i]=parse(Int,x[6]) != 0
end
function read_f64(path,n)
    v=Vector{Float64}(undef,n); open(path,"r") do io; read!(io,v); end; v
end
# Python writes (time, reach) in C order. Julia column-major reshape to
# (reach,time) therefore produces the desired orientation directly.
raw=read_f64(joinpath(work,"qlat.f64"),nreach*ntime)
qlat=reshape(raw,nreach,ntime)
net=RiverNetwork(rid,down)
p=ReachParameters(nreach; length=len,slope=slope,bottom_width=width,side_slope=0.0,
    floodplain_slope=num("floodplain_slope"),bankfull_depth=num("bankfull_depth"),mann_n=num("mann_n"),
    irf_velocity=num("irf_velocity"),irf_diffusivity=num("irf_diffusivity"))
methods = [("irf",IRF()), ("kwt",LagrangianKWT(max_packets=20))]
for (name,method) in methods
    @info "Cameo routing" method=name nreach ntime
    model=RoutingModel(method,net,p;dt=dt,active_reaches=active,headwater_drain_point=2)
    q=route_series(model,qlat)
    open(joinpath(work,"julia_"*name*".f64"),"w") do io; write(io,vec(q)); end
end

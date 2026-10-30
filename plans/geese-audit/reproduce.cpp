#include "barry/barry.hpp"
#include "barry/models/geese.hpp"
#include <iomanip>
using namespace geese;
Geese tree(std::vector<std::vector<size_t>> a) {
 std::vector<size_t> id={0,1,2,3,4}; std::vector<int> pa={-1,0,0,1,1}; std::vector<bool> d(5,true);
 Geese g(a,id,pa,d); counter_gains(g.get_counters(),{0}); counter_loss(g.get_counters(),{0}); g.init(0); return g;
}
int main(int argc,char**argv){
 std::cout<<std::setprecision(16);
 std::vector<double> p={-2,-1,0.3};
 std::string mode=argc>1?argv[1]:"basic";
 if(mode=="counters") {
 int tested=0,failed=0;
 for(bool dup:{false,true})for(size_t et:{0u,1u,2u})for(int parent=0;parent<4;parent++)for(int bits=0;bits<16;bits++){
 PhyloArray a(2,2);a.set_data(new NodeData({1,1},{bool(parent&1),bool(parent&2)},dup),true);
 for(int j=0;j<2;j++)for(int k=0;k<2;k++)a(k,j)=(bits>>(j*2+k))&1;
 PhyloStatsCounter sc(&a);auto cs=sc.get_counters();counter_gains(cs,{0,1},et);counter_loss(cs,{0,1},et);counter_subfun(cs,0,1,et);counter_neofun(cs,0,1,et);counter_neofun_a2b(cs,0,1,et);counter_cogain(cs,0,1,et);
 auto v=sc.count_all();std::vector<double>e(8,0);bool active=et==2||(et==1&&dup)||(et==0&&!dup);
 if(active){for(int k=0;k<2;k++)for(int j=0;j<2;j++){if(!(parent&(1<<k)))e[k]+=a(k,j);else e[k+2]+=1-a(k,j);}
 bool split=(static_cast<size_t>(a(0,0,false))==1&&static_cast<size_t>(a(1,0,false))==0&&static_cast<size_t>(a(0,1,false))==0&&static_cast<size_t>(a(1,1,false))==1)||(static_cast<size_t>(a(0,1,false))==1&&static_cast<size_t>(a(1,1,false))==0&&static_cast<size_t>(a(0,0,false))==0&&static_cast<size_t>(a(1,0,false))==1);
 e[4]=parent==3&&split;e[5]=(parent==1||parent==2)&&split;e[6]=parent==1&&split;if(parent==0)for(int j=0;j<2;j++)e[7]+=a(0,j)*a(1,j);}
 tested++;if(v!=e){failed++;if(failed<4){std::cout<<"counter mismatch parent="<<parent<<" bits="<<bits<<"\n";}}
 }
 std::cout<<"counter configurations="<<tested<<" mismatches="<<failed<<"\n";
 }
 if(mode=="loo-check") {auto g=tree({{9},{9},{1},{0},{1}});auto r=g.predict(p,nullptr,true);for(size_t tip:{2u,3u,4u}){std::vector<std::vector<size_t>>a={{9},{9},{1},{0},{1}};a[tip]={9};auto f=tree(a);auto q=f.predict(p);std::cout<<"tip="<<tip<<" loo="<<r[tip][0]<<" fresh-mask="<<q[tip][0]<<"\n";}}
 if(mode=="cache") {
 auto g=tree({{9},{9},{9},{9},{9}});std::vector<double> a={-8,-8,0.3},b={8,8,0.3};
 for(int i=0;i<100;i++)g.simulate(a);
 g.likelihood(b);g.set_seed(123); int ones=0;for(int i=0;i<10000;i++)ones+=g.simulate(a)[2][0];
 auto f=tree({{9},{9},{9},{9},{9}});f.set_seed(123);int base=0;for(int i=0;i<10000;i++)base+=f.simulate(a)[2][0];
 std::cout<<"simulate(a) after likelihood(b): "<<ones<<"/10000; fresh="<<base<<"/10000\n";
 }
 if(mode=="underflow") {
 size_t n=1100;std::vector<std::vector<size_t>> a(2*n-1,{9});std::vector<size_t> id(2*n-1);std::vector<int> pa(2*n-1,-1);std::vector<bool>d(2*n-1,true);
 for(size_t i=0;i<id.size();++i){id[i]=i;if(i)pa[i]=(i-1)/2;if(i>=n-1)a[i]={1};}
 Geese g(a,id,pa,d);counter_gains(g.get_counters(),{0});counter_loss(g.get_counters(),{0});g.init(0);std::cout<<"1100 fair-coin leaves log="<<g.likelihood({0,0,0},true)<<" expected="<<-1100*std::log(2.)<<"\n";
 }
 if(mode=="hash") {PhyloArray a(17,2), b(17,2);std::vector<bool>x(17,false),y=x;y[0]=true;x[16]=true;y[16]=true;a.set_data(new NodeData({1,1},x),true);b.set_data(new NodeData({1,1},y),true);std::cout<<"distinct parent state keys equal="<<(keygen_full(a,nullptr)==keygen_full(b,nullptr))<<"\n";}
 if(mode=="basic"){
 auto u=tree({{9},{9},{9},{9},{9}});
 std::cout<<"empty probability="<<u.likelihood(p)<<" log="<<u.likelihood(p,true)<<" full-log="<<u.likelihood(p,true,false)<<"\n";
 u.update_annotations(3,{1}); auto f=tree({{9},{9},{9},{1},{9}});
 std::cout<<"new annotation reduced="<<u.likelihood(p)<<" full="<<u.likelihood(p,false,false)<<" fresh="<<f.likelihood(p)<<"\n";
 auto g=tree({{9},{9},{1},{0},{1}});
 auto q=p; q.back()=1000;
 std::cout<<"root logit=1000 likelihood="<<g.likelihood(q)<<"\n";
 auto c=g; c.update_annotations(3,{1});
 auto cf=tree({{9},{9},{1},{1},{1}});
 std::cout<<"copy parent points into source="<<(c.nodes[3].parent==&g.nodes[1])<<"\n";
 std::cout<<"copy update="<<c.likelihood(p,false,false)<<" fresh="<<cf.likelihood(p,false,false)<<"\n";try {std::cout<<"source after copy update="<<g.likelihood(p,false,false)<<"\n";}catch(const std::exception&e){std::cout<<"source error: "<<e.what()<<"\n";}
 }
 if(mode=="internal") {auto g=tree({{1},{9},{1},{0},{1}});std::cout<<"root observed=1 pruning="<<g.likelihood(p)<<" exhaust="<<g.likelihood_exhaust(p)<<"\n";}
 if(mode=="indices") {std::vector<std::vector<size_t>>a={{9},{1},{0}};std::vector<size_t>id={10,20,30};std::vector<int>pa={-1,10,10};std::vector<bool>d(3,true);Geese g(a,id,pa,d);counter_gains(g.get_counters(),{0});counter_loss(g.get_counters(),{0});g.init(0);std::cout<<"pruning="<<g.likelihood(p)<<"\n";try{std::cout<<"exhaust="<<g.likelihood_exhaust(p)<<"\n";}catch(const std::exception&e){std::cout<<"exhaust error="<<e.what()<<"\n";}}
 if(mode=="loo") {auto g=tree({{9},{9},{1},{0},{1}}); auto r=g.predict(p,nullptr,true);for(auto v:r)std::cout<<v[0]<<" ";std::cout<<"\n";}
 if(mode=="update") { auto g=tree({{9},{9},{1},{0},{1}});g.update_annotations(3,{9});g.update_annotations(4,{9});auto f=tree({{9},{9},{1},{9},{9}});std::cout<<"updated="<<g.likelihood(p)<<" fresh="<<f.likelihood(p)<<"\n"; }
 if(mode=="singleton") {std::vector<std::vector<size_t>> a={{1}};std::vector<size_t> id={0};std::vector<int> pa={-1};std::vector<bool>d={true};Geese g(a,id,pa,d);counter_gains(g.get_counters(),{0});g.init(0);std::cout<<g.likelihood({-2,0.3})<<"\n";}
}

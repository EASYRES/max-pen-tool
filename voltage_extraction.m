function [Voltage]=voltage_extraction(BusNames,NodeNames,k)


lb=length(BusNames);
ln=length(NodeNames);

for x=1:lb
    node1=[BusNames{x},'.1'];
    node2=[BusNames{x},'.2'];
    node3=[BusNames{x},'.3'];
    node4=[BusNames{x},'.4'];
    VNn(x)=0;
   
    for j=1:ln
        pha=strcmp(node1,NodeNames{j});
        if pha==1
            VAn(x)=k((2*j)-1)+1i*k(2*j);
        end
    end
    
    for j=1:ln
        phb=strcmp(node2,NodeNames{j});
        if phb==1
            VBn(x)=k((2*j)-1)+1i*k(2*j);
        end
    end
    
    for j=1:ln
        phc=strcmp(node3,NodeNames{j});
        if phc==1
            VCn(x)=k((2*j)-1)+1i*k(2*j);
        end
    end
    
    for j=1:ln
        neutral=strcmp(node4,NodeNames{j});
        if neutral==1
            VNn(x)=k((2*j)-1)+1i*k(2*j);
        end
    end
    
end

VAN_vec=VAn-VNn;
VBN_vec=VBn-VNn;
VCN_vec=VCn-VNn;

Voltage.VAN_vec=VAN_vec;
Voltage.VBN_vec=VBN_vec;
Voltage.VCN_vec=VCN_vec;

Voltage.V1=sqrt(3)*abs((VAN_vec+VBN_vec*exp(1i*2*pi/3)+VCN_vec*exp(1i*4*pi/3))/3);
Voltage.V1_mag_pu(1)=Voltage.V1(1)/(20000);
Voltage.V1_mag_pu=[Voltage.V1_mag_pu(1),Voltage.V1(2:end)/(400)];
Voltage.Vmax=sqrt(3)*max([abs(VAN_vec);abs(VBN_vec);abs(VCN_vec)]);
Voltage.Vmax_pu(1)=Voltage.Vmax(1)/(20000);
Voltage.Vmax_pu=[Voltage.Vmax_pu(1),Voltage.Vmax(2:end)/400];

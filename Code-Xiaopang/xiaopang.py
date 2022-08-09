#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import math
import numpy as np

class Quaternary_huffman:
    def __init__(self,**info):
        self.data=info['data']

    def codegenerate(self):
        quat_huf_init=self.Quaternary_huffman_init(data=self.data)
        symbol_composition, symbol_frequency=quat_huf_init.p_generate()
        seat,codewords=quat_huf_init.seat_and_codewords(symbol_frequency)
        
        code=''
        for data_current in self.data:
            index=symbol_composition.index(data_current)
            code=code+codewords[index]
        return code,symbol_composition,seat
        
    class Quaternary_huffman_init():
        def __init__(self,**info):
            self.data=info['data']
        
        def p_generate(self):
            
            symbol_composition=[]
            symbol_frequency=[]     
            
            for data_current in self.data:
                t_existence = data_current in symbol_composition
                if t_existence != True:
                    symbol_composition.append(data_current)     
                    
            for symbol_current in symbol_composition:
                symbol_frequency.append(self.data.count(symbol_current))
                
            record=[symbol_composition,symbol_frequency]
            record=list(record)
            record=list(map(list, zip(*record)))	
            
            def takeSecond(elem):
                return elem[1]
            record.sort(key=takeSecond)
            
            symbol_composition=[]
            symbol_frequency=[]
             
            for c in record:
                symbol_composition.append(c[0])
                symbol_frequency.append(c[1]/len(self.data))
                
            return symbol_composition, symbol_frequency
        
        def seat_and_codewords(self,symbol_frequency):
            n=(len(symbol_frequency)-1)%3;
            if(n==0):
                symbol_frequency_fictitious=symbol_frequency;
            elif(n==1):
                symbol_frequency_fictitious=[0,0]+symbol_frequency;
            elif(n==2):
                symbol_frequency_fictitious=[0]+symbol_frequency;
            
            symbol_frequency_fictitious_initial=symbol_frequency_fictitious
            seat=[]
            while(len(symbol_frequency_fictitious)!=4): 
                smallest=sum(symbol_frequency_fictitious[0:4]) 
                p_sort_before=np.append(smallest,symbol_frequency_fictitious[4:])  
                p_sort_after=np.sort(p_sort_before)     
                seat.append(np.where(p_sort_after==smallest)[0][0]);    
                symbol_frequency_fictitious=p_sort_after
        
            table_quaternary_value=['1','2','3','4']
            codewords=[];[codewords.append(quaternary_value) for quaternary_value in table_quaternary_value]
            
            codewords_temp=codewords
            for i in range(1,len(seat)+1):
                seat_c=seat[-i]
                codewords=[];[codewords.append(codewords_temp[seat_c]+quaternary_value) for quaternary_value in table_quaternary_value]
                codewords_temp.pop(seat_c)
                codewords=codewords+codewords_temp
                codewords_temp=codewords
                
            codewords=codewords[0+len(symbol_frequency_fictitious_initial)-len(symbol_frequency):]    
        
            return seat,codewords

class Code_table_addition:
    def __init__(self,**info):
        self.symbol_composition=info['symbol_composition']
        self.seat=info['seat']

    def f_huffman_table(self):
        information_table=self.f_information_table()
        tree_table=self.f_tree_table()
        header1=self.f_header1()
        header2=self.f_header2()
        huffman_table=header1+header2+information_table+tree_table
        return huffman_table
    
    def f_information_table(self):        
        information_table=[]
        for word in self.symbol_composition:
            information_table.append(ord(word))     
        information_table=self.f_decimal_to_quaternary(information_table)
        return information_table
    
    def f_tree_table(self):
        tree_table=self.f_decimal_to_quaternary(self.seat)
        return tree_table
    
    def f_header1(self):
        header1=self.f_decimal_to_quaternary([len(self.symbol_composition)])        
        return header1

    def f_header2(self):
        header2=self.f_decimal_to_quaternary([len(self.seat)])        
        return header2

    def f_decimal_to_quaternary(self,data):
        rec_quaternary=[]      
        for current in data:
            value_2=str(bin(current)[2:].zfill(8));       
            value_4=''             
            for i in range(0,8,2):  
                if value_2[i]=='0' and  value_2[i+1]=='0': value_4=value_4+'1' 
                if value_2[i]=='0' and  value_2[i+1]=='1': value_4=value_4+'2' 
                if value_2[i]=='1' and  value_2[i+1]=='0': value_4=value_4+'3' 
                if value_2[i]=='1' and  value_2[i+1]=='1': value_4=value_4+'4'         
            rec_quaternary.append(value_4)           
        rec_quaternary=''.join(rec_quaternary)   
        return rec_quaternary
     
class Extension_mapping:
    def __init__(self,**info):
        [self.code5,self.code6]=c5c6()
        self.code_initial=info['code_initial']
        
    def f_extension_mapping_twice(self):
        extension_first=self.code_5to6(self.code_initial)
        extension_second=self.encode_homopolymer(extension_first)
        return extension_second

    def code_5to6(self,code):
        def anyToDecimal(num,n):
           baseStr = {"0":0,"1":1,"2":2,"3":3,"4":4,"5":5,"6":6,"7":7,"8":8,"9":9,
                      "a":10,"b":11,"c":12,"d":13,"e":14,"f":15,"g":16,"h":17,"i":18,"j":19}
           new_num = 0
           nNum = len(num) - 1
           for i in num:
               new_num = new_num  + baseStr[i]*pow(n,nNum)
               nNum = nNum -1
           return new_num
       
        more=len(code)%5
        temp='1234'
        rec=''
        for i in range(math.floor(len(code)/5)):
            mid=str(int(code[i*5:(i+1)*5])-11111)
            locc=anyToDecimal(mid,4)
            rec=rec+str(self.code6[locc])
        
        if more!=0:
            rec=rec+code[-more:]+temp[0:(5-more)]
        return rec
        
    def encode_homopolymer(self,code):
        rec=[]
        i=0
        while(i<=len(code)-1):
            if code[i:i+4]=='4444':
                rec.append('12344');i=i+4
            
            elif code[i:i+4]=='3333':
                rec.append('12343');i=i+4
            
            elif code[i:i+4]=='2222':
                rec.append('12342');i=i+4
            
            elif code[i:i+4]=='1111':
                rec.append('12341');i=i+4
            
            else:
                rec.append(code[i]);i=i+1
            
        code=''.join(rec)
        return code

class Segmentation_index(Extension_mapping):
    def __init__(self,**info):
        [self.code5,self.code6]=c5c6()
        self.sequences_initial=info['sequences_initial']
        self.n_bases_payload=info['n_bases_payload']
        
    def num_to_N(self,sequences_silicon):
        sequences=[]
        for dna in sequences_silicon:
            rec=[]
            for N in dna:
                if N=='1':rec.append('A')
                if N=='2':rec.append('T')
                if N=='3':rec.append('G')
                if N=='4':rec.append('C')
            sequences.append(''.join(rec))
        return sequences
    
    def dna_generate(self):
        
        num_dna=math.ceil(len(self.sequences_initial)/self.n_bases_payload)    
        
        rec_indices=[]; 
        rec_payload=[]; 
        
        if num_dna>4**5: 
            value_threshold=20
        else:
            value_threshold=10
            
        for ic in range(num_dna):
            mid1=bin(ic);mid2=mid1[2:]
            while len(mid2)!=value_threshold:
                mid2='0'+mid2
        
            mid_need=''
            for i in range(0,10,2):
                if mid2[i]=='0' and  mid2[i+1]=='0':
                    mid_need=mid_need+'1'
                if mid2[i]=='0' and  mid2[i+1]=='1':
                    mid_need=mid_need+'2'
                if mid2[i]=='1' and  mid2[i+1]=='0':
                    mid_need=mid_need+'3'
                if mid2[i]=='1' and  mid2[i+1]=='1':
                    mid_need=mid_need+'4'
    
            mid_need=self.code_5to6(mid_need)    
            
            rec_indices.append(mid_need)
            rec_payload.append(self.sequences_initial[ic*self.n_bases_payload:(ic+1)*self.n_bases_payload])
    
        rec_sequence_number=[];
        for i in range(num_dna):
            rec_sequence_number.append(rec_indices[i]+rec_payload[i])
        
        return rec_sequence_number
    
    def f_index_and_num2base(self):       
        sequences_silicon=self.dna_generate()      
        dna_sequences=self.num_to_N(sequences_silicon)
        
        return dna_sequences
     
def c5c6():
    code6=[322424, 323241, 344121, 314131, 233213, 113342, 421324, 131244, 142313, 133421, 143132, 241441, 321243, 142431, 243232, 324124, 431131, 311342, 123331, 231324, 324231, 424132, 143242, 312332, 121443, 411341, 131323, 314321, 144131, 124133, 232342, 341311, 113233, 233311, 132332, 311331, 312431, 421342, 231423, 231313, 134311, 413241, 214233, 132442, 311424, 322442, 131331, 421143, 342231, 131414, 131342, 134114, 214324, 342113, 241323, 231332, 341131, 144223, 331411, 321144, 414214, 331141, 214341, 223431, 313113, 421431, 422441, 424321, 323224, 231442, 141432, 143422, 342141, 232233, 123144, 243124, 342242, 314241, 321441, 314224, 323114, 342321, 424141, 144411, 332241, 224423, 312244, 241243, 221433, 122433, 232244, 412424, 422323, 324242, 132324, 312324, 231144, 242331, 214244, 342124, 344211, 241414, 114243, 143321, 312414, 244224, 423131, 134224, 431214, 431241, 133232, 431113, 423242, 422234, 223332, 224314, 413421, 424113, 113134, 242314, 213314, 213133, 213424, 133223, 341422, 244114, 414311, 212433, 313141, 313231, 213331, 414231, 312423, 242234, 312134, 313214, 412313, 244213, 233242, 314411, 323142, 141314, 243131, 413114, 414242, 332224, 331124, 331242, 214442, 342223, 244231, 142144, 321423, 233422, 313422, 341241, 423214, 123243, 144113, 231431, 411314, 323232, 341214, 133214, 432132, 421414, 243241, 213441, 211344, 243322, 243411, 213432, 143231, 231133, 124431, 421244, 223134, 232414, 244141, 223442, 144142, 244321, 322133, 131441, 324422, 321432, 214313, 142342, 422314, 341142, 144232, 411243, 411423, 331322, 223323, 114314, 332142, 113413, 133131, 143224, 322313, 242413, 214134, 421134, 414422, 431421, 113323, 413224, 114432, 124234, 112443, 143114, 331213, 331132, 144214, 422332, 131233, 133322, 322341, 133113, 142233, 232441, 223143, 324223, 314232, 142414, 124442, 232432, 422413, 121334, 233141, 144124, 313223, 112334, 224432, 132143, 314213, 213144, 223413, 223341, 113431, 124323, 224243, 242424, 311413, 323411, 143141, 412143, 313132, 132341, 214143, 141134, 224144, 232323, 411432, 332311, 223342, 423142, 133321, 213233, 123323, 213332, 142424, 242341, 124342, 233142, 314322, 214234, 244113, 241431, 422244, 422341, 341232, 243132, 422143, 241144, 131134, 214323, 132234, 421314, 221334, 144132, 421413, 242233, 231342, 213442, 321143, 114341, 232313, 244422, 223234, 123143, 211434, 114442, 411413, 113243, 324321, 322332, 314141, 311432, 213143, 221443, 322144, 223243, 323311, 123341, 314231, 243223, 424322, 332114, 423213, 421432, 332131, 424241, 132243, 321313, 414321, 242432, 232424, 233322, 422442, 424311, 314223, 314114, 431242, 432141, 414142, 121344, 131341, 341213, 212334, 422324, 341321, 133132, 124424, 112434, 324411, 322414, 242423, 214331, 311323, 331113, 321342, 124441, 144213, 144322, 422233, 312413, 342132, 412441, 421133, 413242, 424232, 243242, 121433, 421441, 423411, 412323, 321244, 114431, 324113, 132432, 431132, 413141, 412432, 422423, 122343, 124432, 223331, 332213, 114324, 411324, 411332, 413223, 133422, 324142, 224341, 421323, 141313, 344122, 223423, 312143, 314132, 231233, 341411, 142442, 142314, 214133, 342421, 213313, 142244, 414124, 143113, 423231, 331142, 133124, 412331, 132133, 423224, 123244, 331311, 331214, 323141, 313142, 412133, 134214, 412414, 122334, 233224, 142332, 342241, 243214, 332321, 142323, 413214, 224244, 421144, 341224, 233114, 342322, 413322, 334221, 133241, 323132, 422431, 224143, 321233, 324214, 213244, 113331, 311243, 231414, 141441, 144311, 413132, 242143, 323322, 114313, 422313, 331231, 233232, 124331, 431322, 114244, 114423, 324131, 133411, 224134, 244124, 413311, 411233, 231143, 313321, 214424, 141414, 331422, 424421, 332242, 144241, 241133, 431124, 411313, 312314, 123134, 223414, 224442, 224324, 314242, 313411, 224313, 423114, 421331, 313124, 134142, 313131, 341114, 423422, 114134, 132423, 214342, 343221, 232332, 134421, 132331, 241413, 134231, 133213, 113234, 241314, 141143, 242324, 311441, 123314, 124144, 432113, 223133, 314421, 213324, 412234, 124233, 422134, 342214, 412144, 431223, 414113, 332411, 321134, 411342, 421424, 232134, 313213, 243421, 134322, 141342, 424422, 143142, 321414, 214314, 143223, 422243, 312243, 244142, 214431, 423232, 224234, 213234, 413321, 114323, 232442, 311244, 414131, 323213, 124243, 412332, 142341, 322143, 431321, 313421, 311341, 141233, 242342, 241341, 124332, 432131, 322342, 313114, 312441, 133311, 324132, 224414, 113341, 144141, 423113, 314214, 232133, 423311, 423124, 241313, 131424, 133141, 243141, 312313, 123233, 412431, 143214, 312432, 432142, 232243, 231323, 244232, 131314, 144421, 342311, 414232, 132313, 113133, 244223, 321324, 243213, 211443, 114441, 421233, 241424, 224323, 331421, 213423, 313241, 341223, 242441, 422331, 411441, 424242, 431422, 322323, 214332, 214144, 412324, 322244, 233132, 411414, 144224, 133224, 411331, 132134, 132414, 144242, 141423, 134232, 323131, 213243, 232413, 131133, 132244, 114233, 221343, 343211, 342114, 133242, 224342, 323242, 421341, 213414, 431141, 134213, 121434, 141244, 244421, 322314, 141324, 431114, 332132, 242323, 314311, 131332, 422432, 114143, 311143, 223244, 123324, 332223, 413213, 244131, 424231, 131413, 331224, 414322, 233223, 321133, 244241, 223424, 123313, 432213, 424213, 233124, 311431, 412244, 314142, 323321, 231441, 231134, 231331, 134124, 332231, 334121, 421332, 113424, 233421, 311314, 421313, 412233, 242134, 313232, 332214, 314124, 241332, 421423, 241442, 233113, 413231, 322431, 113442, 143411, 412423, 424224, 423141, 141331, 124413, 142133, 113144, 322134, 231432, 334211, 144321, 223314, 243321, 412442, 131442, 432124, 144422, 134131, 232234, 331232, 322413, 244311, 223432, 423421, 422424, 123342, 232341, 231244, 114342, 423223, 312234, 313224, 414241, 431232, 424114, 114414, 321331, 142234, 224133, 332422, 311442, 122443, 311233, 324421, 324241, 134241, 232431, 134411, 421442, 343121, 322233, 342232, 233231, 332141, 431411, 131144, 244322, 232144, 132441, 312331, 213342, 341132, 311332, 131243, 314113, 242244, 131431, 312424, 124341, 142243, 414223, 113314, 224331, 312133, 411134, 124314, 423322, 314422, 112433, 323422, 114331, 113332, 214413, 414411, 232314, 142413, 142432, 143213, 124244, 313322, 142423, 431311, 131143, 213413, 244411, 431224, 221434, 124313, 243422, 132413, 421243, 141424, 312442, 334122, 223441, 422414, 414141, 143232, 224233, 141442, 123133, 242144, 321424, 131432, 141144, 242243, 412243, 114413, 224424, 324213, 321323, 321442, 423132, 224332, 422342, 124423, 413142, 412134, 232331, 242332, 231424, 311414, 312323, 241342, 324141, 321413, 323231, 114424, 224431, 241134, 332421, 232423, 142143, 221344, 413113, 113313, 142134, 324322, 311324, 422133, 321431, 131313, 124143, 223144, 114133, 313242, 431231, 134321, 432114, 211343, 323124, 214243, 213341, 241233, 311423, 224413, 231243, 114332, 414224, 124324, 141243, 211433, 342422, 124414, 213323, 113244, 411323, 243113, 134242, 324224, 132431, 244132, 143311, 342131, 341141, 424131, 332232, 231314, 342411, 142324, 423321, 241324, 243142, 242442, 134113, 131324, 243224, 223313, 331114, 211334, 341242, 233241, 322331, 241423, 324232, 321332, 241244, 143131, 242414, 413422, 141332, 213431, 242313, 324311, 232324, 331223, 323113, 223324, 144114, 412314, 132144, 132323, 233321, 212434, 414421, 231341, 232143, 113441, 324114, 142441, 342224, 143124, 122434, 113414, 123234, 132314, 323223, 144231, 332124, 321341, 343122, 132233, 411244, 413131, 113432, 313311, 322441, 141413, 321314, 423241, 312144, 411424, 431213, 113423, 413232, 143241, 322324, 424411, 322432, 341231, 331131, 344221, 114234, 242133, 331241, 413411, 341113, 214441, 322234, 414213, 341124, 233214, 312233, 233131, 143322, 431142, 341421, 133231, 223233, 123332, 341322, 331321, 244214, 214432, 134223, 141341, 133114, 411143, 134132, 141323, 122344, 412413, 214423, 224441, 411431, 424214, 114144, 342142, 143421, 121343, 332322, 134422, 131423, 212443, 241143, 141431, 132342, 424124, 323214, 323421, 142331, 422144, 124134, 241331, 134141, 132424, 322243, 332113, 141133, 424223, 413124, 414132, 233411, 243231, 214414, 322423, 113324, 411442, 242431, 244242, 243311, 231413, 133142, 213134, 424142, 243114, 311313, 241432, 342213, 414114, 113143, 311134]            
    code5=[11111, 11112, 11113, 11114, 11121, 11122, 11123, 11124, 11131, 11132, 11133, 11134, 11141, 11142, 11143, 11144, 11211, 11212, 11213, 11214, 11221, 11222, 11223, 11224, 11231, 11232, 11233, 11234, 11241, 11242, 11243, 11244, 11311, 11312, 11313, 11314, 11321, 11322, 11323, 11324, 11331, 11332, 11333, 11334, 11341, 11342, 11343, 11344, 11411, 11412, 11413, 11414, 11421, 11422, 11423, 11424, 11431, 11432, 11433, 11434, 11441, 11442, 11443, 11444, 12111, 12112, 12113, 12114, 12121, 12122, 12123, 12124, 12131, 12132, 12133, 12134, 12141, 12142, 12143, 12144, 12211, 12212, 12213, 12214, 12221, 12222, 12223, 12224, 12231, 12232, 12233, 12234, 12241, 12242, 12243, 12244, 12311, 12312, 12313, 12314, 12321, 12322, 12323, 12324, 12331, 12332, 12333, 12334, 12341, 12342, 12343, 12344, 12411, 12412, 12413, 12414, 12421, 12422, 12423, 12424, 12431, 12432, 12433, 12434, 12441, 12442, 12443, 12444, 13111, 13112, 13113, 13114, 13121, 13122, 13123, 13124, 13131, 13132, 13133, 13134, 13141, 13142, 13143, 13144, 13211, 13212, 13213, 13214, 13221, 13222, 13223, 13224, 13231, 13232, 13233, 13234, 13241, 13242, 13243, 13244, 13311, 13312, 13313, 13314, 13321, 13322, 13323, 13324, 13331, 13332, 13333, 13334, 13341, 13342, 13343, 13344, 13411, 13412, 13413, 13414, 13421, 13422, 13423, 13424, 13431, 13432, 13433, 13434, 13441, 13442, 13443, 13444, 14111, 14112, 14113, 14114, 14121, 14122, 14123, 14124, 14131, 14132, 14133, 14134, 14141, 14142, 14143, 14144, 14211, 14212, 14213, 14214, 14221, 14222, 14223, 14224, 14231, 14232, 14233, 14234, 14241, 14242, 14243, 14244, 14311, 14312, 14313, 14314, 14321, 14322, 14323, 14324, 14331, 14332, 14333, 14334, 14341, 14342, 14343, 14344, 14411, 14412, 14413, 14414, 14421, 14422, 14423, 14424, 14431, 14432, 14433, 14434, 14441, 14442, 14443, 14444, 21111, 21112, 21113, 21114, 21121, 21122, 21123, 21124, 21131, 21132, 21133, 21134, 21141, 21142, 21143, 21144, 21211, 21212, 21213, 21214, 21221, 21222, 21223, 21224, 21231, 21232, 21233, 21234, 21241, 21242, 21243, 21244, 21311, 21312, 21313, 21314, 21321, 21322, 21323, 21324, 21331, 21332, 21333, 21334, 21341, 21342, 21343, 21344, 21411, 21412, 21413, 21414, 21421, 21422, 21423, 21424, 21431, 21432, 21433, 21434, 21441, 21442, 21443, 21444, 22111, 22112, 22113, 22114, 22121, 22122, 22123, 22124, 22131, 22132, 22133, 22134, 22141, 22142, 22143, 22144, 22211, 22212, 22213, 22214, 22221, 22222, 22223, 22224, 22231, 22232, 22233, 22234, 22241, 22242, 22243, 22244, 22311, 22312, 22313, 22314, 22321, 22322, 22323, 22324, 22331, 22332, 22333, 22334, 22341, 22342, 22343, 22344, 22411, 22412, 22413, 22414, 22421, 22422, 22423, 22424, 22431, 22432, 22433, 22434, 22441, 22442, 22443, 22444, 23111, 23112, 23113, 23114, 23121, 23122, 23123, 23124, 23131, 23132, 23133, 23134, 23141, 23142, 23143, 23144, 23211, 23212, 23213, 23214, 23221, 23222, 23223, 23224, 23231, 23232, 23233, 23234, 23241, 23242, 23243, 23244, 23311, 23312, 23313, 23314, 23321, 23322, 23323, 23324, 23331, 23332, 23333, 23334, 23341, 23342, 23343, 23344, 23411, 23412, 23413, 23414, 23421, 23422, 23423, 23424, 23431, 23432, 23433, 23434, 23441, 23442, 23443, 23444, 24111, 24112, 24113, 24114, 24121, 24122, 24123, 24124, 24131, 24132, 24133, 24134, 24141, 24142, 24143, 24144, 24211, 24212, 24213, 24214, 24221, 24222, 24223, 24224, 24231, 24232, 24233, 24234, 24241, 24242, 24243, 24244, 24311, 24312, 24313, 24314, 24321, 24322, 24323, 24324, 24331, 24332, 24333, 24334, 24341, 24342, 24343, 24344, 24411, 24412, 24413, 24414, 24421, 24422, 24423, 24424, 24431, 24432, 24433, 24434, 24441, 24442, 24443, 24444, 31111, 31112, 31113, 31114, 31121, 31122, 31123, 31124, 31131, 31132, 31133, 31134, 31141, 31142, 31143, 31144, 31211, 31212, 31213, 31214, 31221, 31222, 31223, 31224, 31231, 31232, 31233, 31234, 31241, 31242, 31243, 31244, 31311, 31312, 31313, 31314, 31321, 31322, 31323, 31324, 31331, 31332, 31333, 31334, 31341, 31342, 31343, 31344, 31411, 31412, 31413, 31414, 31421, 31422, 31423, 31424, 31431, 31432, 31433, 31434, 31441, 31442, 31443, 31444, 32111, 32112, 32113, 32114, 32121, 32122, 32123, 32124, 32131, 32132, 32133, 32134, 32141, 32142, 32143, 32144, 32211, 32212, 32213, 32214, 32221, 32222, 32223, 32224, 32231, 32232, 32233, 32234, 32241, 32242, 32243, 32244, 32311, 32312, 32313, 32314, 32321, 32322, 32323, 32324, 32331, 32332, 32333, 32334, 32341, 32342, 32343, 32344, 32411, 32412, 32413, 32414, 32421, 32422, 32423, 32424, 32431, 32432, 32433, 32434, 32441, 32442, 32443, 32444, 33111, 33112, 33113, 33114, 33121, 33122, 33123, 33124, 33131, 33132, 33133, 33134, 33141, 33142, 33143, 33144, 33211, 33212, 33213, 33214, 33221, 33222, 33223, 33224, 33231, 33232, 33233, 33234, 33241, 33242, 33243, 33244, 33311, 33312, 33313, 33314, 33321, 33322, 33323, 33324, 33331, 33332, 33333, 33334, 33341, 33342, 33343, 33344, 33411, 33412, 33413, 33414, 33421, 33422, 33423, 33424, 33431, 33432, 33433, 33434, 33441, 33442, 33443, 33444, 34111, 34112, 34113, 34114, 34121, 34122, 34123, 34124, 34131, 34132, 34133, 34134, 34141, 34142, 34143, 34144, 34211, 34212, 34213, 34214, 34221, 34222, 34223, 34224, 34231, 34232, 34233, 34234, 34241, 34242, 34243, 34244, 34311, 34312, 34313, 34314, 34321, 34322, 34323, 34324, 34331, 34332, 34333, 34334, 34341, 34342, 34343, 34344, 34411, 34412, 34413, 34414, 34421, 34422, 34423, 34424, 34431, 34432, 34433, 34434, 34441, 34442, 34443, 34444, 41111, 41112, 41113, 41114, 41121, 41122, 41123, 41124, 41131, 41132, 41133, 41134, 41141, 41142, 41143, 41144, 41211, 41212, 41213, 41214, 41221, 41222, 41223, 41224, 41231, 41232, 41233, 41234, 41241, 41242, 41243, 41244, 41311, 41312, 41313, 41314, 41321, 41322, 41323, 41324, 41331, 41332, 41333, 41334, 41341, 41342, 41343, 41344, 41411, 41412, 41413, 41414, 41421, 41422, 41423, 41424, 41431, 41432, 41433, 41434, 41441, 41442, 41443, 41444, 42111, 42112, 42113, 42114, 42121, 42122, 42123, 42124, 42131, 42132, 42133, 42134, 42141, 42142, 42143, 42144, 42211, 42212, 42213, 42214, 42221, 42222, 42223, 42224, 42231, 42232, 42233, 42234, 42241, 42242, 42243, 42244, 42311, 42312, 42313, 42314, 42321, 42322, 42323, 42324, 42331, 42332, 42333, 42334, 42341, 42342, 42343, 42344, 42411, 42412, 42413, 42414, 42421, 42422, 42423, 42424, 42431, 42432, 42433, 42434, 42441, 42442, 42443, 42444, 43111, 43112, 43113, 43114, 43121, 43122, 43123, 43124, 43131, 43132, 43133, 43134, 43141, 43142, 43143, 43144, 43211, 43212, 43213, 43214, 43221, 43222, 43223, 43224, 43231, 43232, 43233, 43234, 43241, 43242, 43243, 43244, 43311, 43312, 43313, 43314, 43321, 43322, 43323, 43324, 43331, 43332, 43333, 43334, 43341, 43342, 43343, 43344, 43411, 43412, 43413, 43414, 43421, 43422, 43423, 43424, 43431, 43432, 43433, 43434, 43441, 43442, 43443, 43444, 44111, 44112, 44113, 44114, 44121, 44122, 44123, 44124, 44131, 44132, 44133, 44134, 44141, 44142, 44143, 44144, 44211, 44212, 44213, 44214, 44221, 44222, 44223, 44224, 44231, 44232, 44233, 44234, 44241, 44242, 44243, 44244, 44311, 44312, 44313, 44314, 44321, 44322, 44323, 44324, 44331, 44332, 44333, 44334, 44341, 44342, 44343, 44344, 44411, 44412, 44413, 44414, 44421, 44422, 44423, 44424, 44431, 44432, 44433, 44434, 44441, 44442, 44443, 44444]    
    return code5,code6  
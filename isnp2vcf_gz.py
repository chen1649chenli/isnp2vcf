#!/pbp/local/bin/python3
import os;
import math;
import gzip;
import shutil;
def convert_nonvariant_line(array):
    '''(array) -> string
       This function convert a filed array that is derived from a nonvariant isnp line into a line of string in vcf format
    '''

    chrom = array[0];
    pos = array[1];
    id = ".";
    ref = array[2];
    alt = ".";
    qual = array[7];
    filter = "PASS";
    info = "DP=" + array[11];
    format = "GT:DP";
    genotype ="0/0:" + array[6]

    vcf_line = chrom + "\t" + pos + "\t" + id + "\t" + ref + "\t" + alt + "\t" + qual + "\t" + filter + "\t" + info + "\t" + format + "\t" + genotype + "\n";
    return vcf_line;



def convert_variant_line(array):
    '''(array) -> string
       This function convert a filed array that is derived from a variant isnp line into a line of string in vcf format
    '''

    chrom = array[0];
    pos = array[1];
    id = ".";
    ref = array[2];
    alt = get_alt(array[3],array[2]);
    qual = array[7];
    filter = "PASS";
    info = "DP=" + array[11];
    format = "GT:DP:GQ";
    genotype = get_genotype(array[3],array[2],array[11],array[4]);


    vcf_line = chrom + "\t" + pos + "\t" + id + "\t" + ref + "\t" + alt + "\t" + qual + "\t" + filter + "\t" + info + "\t" + format + "\t" + genotype + "\n";
    return vcf_line;

def get_alt(alt,ref):
    '''(string) -> string
       This function convert alternative allele in isnp format to vcf format
    '''
    ref = ref.upper();
    alt = alt.upper();
    alt_array = [];
    if alt in ["A","T","C","G"]:
        alt_array.append(alt);
    elif alt ==  "W":
        alt_array = ["A","T"];
    elif alt == "S":
        alt_array = ["C","G"];
    elif alt == "M":
        alt_array = ["A","C"];
    elif alt == "K":
        alt_array = ["G","T"];
    elif alt == "R":
        alt_array = ["A","G"];
    elif alt == "Y":
        alt_array = ["C","T"];
    else:
        alt_array = "N";
    if ref in alt_array:
        alt_array.remove(ref);
    alt = str.join(",",alt_array);

    return alt;

def get_genotype(alt,ref,depth,pValue):

    '''(string) -> string
       This function calculates the genotype value for the genotype column in the vcf files
    '''
    GQ = float(pValue)
    if GQ == 0:
        GQ = 0.0000001;
    GQ = abs((math.log10(GQ))*(-10));
    GQ = str(GQ);
    GQ = GQ[:4];

    DP = depth;

    ref = ref.upper();
    alt = alt.upper();
    alt_array = [];
    if alt in ["A","T","C","G"]:
        GT = "1/1";
    elif alt == "W":
        alt_array = ["A","T"];
        if ref in alt_array:
            GT = "0/1";
        else:
            GT = "1/2";
    elif alt == "S":
        alt_array = ["C","G"];
        if ref in alt_array:
            GT = "0/1";
        else:
            GT = "1/2";
    elif alt == "M":
        alt_array = ["A","C"];
        if ref in alt_array:
            GT = "0/1";
        else:
            GT = "1/2";
    elif alt == "K":
        alt_array = ["G","T"];
        if ref in alt_array:
            GT = "0/1";
                    else:
            GT = "1/2";
    elif alt == "R":
        alt_array = ["A","G"];
        if ref in alt_array:
            GT = "0/1";
        else:
            GT = "1/2";
    elif alt == "Y":
        alt_array = ["C","T"];
        if ref in alt_array:
            GT = "0/1";
        else:
            GT = "1/2";
    else:
        alt_array = "N";
        GT = "0/0"

    genotype = GT + ":" + DP + ":" + GQ;
    return genotype;
    
    vcf_output_dir = "/import/scratch/user/llchen5/isnp2vcf/final_data";
current_dir = os.getcwd();
with open ("isnp_file","r")as isnp_source_file:
    for file_line in isnp_source_file:
        file_array = file_line.split('\t');
        file_path = file_array[0];
        seqLIMS_sample = file_array[1];
        seqLIMS_project = file_array[2];
        inventoryBID = file_array[3];
        sequence_technology = file_array[4];
        capture_probe_set = file_array[5];
        seqoID = "GenRe_"+inventoryBID;
        pedigree = file_array[6];
        header = "##fileformat=VCFv4.1\n"+"##Reference_Genome_Business_Key=<Species=Glycine max,Subspecies=Williams82,Source=Monsanto,Version=V1>\n"+"##SeqLIMS_Sample="+ seqLIMS_sample + "\n##SeqLIMS_Project=" + seqLIMS_project + "\n##Capture_Probe_Set=" + capture_probe_set + "\n##Sequencing_Technology=" + sequence_technology + "\n##InventoryBID=" + inventoryBID + "\n##SeqoID=" + seqoID +"\n#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT\t" + pedigree + "\n";

        tmp_output_dir = current_dir + "/" + seqLIMS_sample.replace("-","_");
        if os.path.isdir(tmp_output_dir):
            print("folder " + tmp_output_dir + " exist, please remove it first\n");
            continue;
        else:
            os.mkdir(tmp_output_dir);
            vcf_file_path = tmp_output_dir +"/" + seqLIMS_sample.replace("-","_") + ".vcf";
            vcf_zipped_file = vcf_file_path + ".gz";
            finished_file = "/import/scratch/user/llchen5/isnp2vcf/final_data/" +seqLIMS_sample.replace("-","_") + ".vcf.gz";
            vcf_file = open (vcf_file_path,"w")
            vcf_file.write(header);

        with gzip.open (file_path) as isnp_file:
            for line in isnp_file:
                line = line.decode('utf-8');
                line =  line.rstrip();
                if line.startswith("#"):
                    continue;
                field_array = line.split("\t");
                if len(field_array[2]) > 1:
                    continue;
                if field_array[2] == field_array[3]:
                   vcf_file.write(convert_nonvariant_line(field_array));
                else:
                   vcf_file.write(convert_variant_line(field_array));

        vcf_file.close();

        with open(vcf_file_path, 'rb') as f_in:
            with gzip.open(vcf_zipped_file, 'wb') as f_out:
                f_out.writelines(f_in);

        os.remove(vcf_file_path);
        shutil.move(vcf_zipped_file,finished_file);






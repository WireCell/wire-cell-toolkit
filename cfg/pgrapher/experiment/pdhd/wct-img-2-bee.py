import sys, os, glob

def main(fp1,fp2,fp3,fp4):
    if (os.path.exists('data/0')):
        print('found old data, removing ...')
        os.system('rm -rf data')
    if (os.path.exists('upload.zip')):
        os.system('rm -f upload.zip')
    os.system('mkdir -p data/0')
    density = 0.1
    cmd = 'wirecell-img bee-blobs -g protodunehd -s uniform -d %f --speed "-1.6*mm/us" --t0 "250*us" --x0 "-  358*cm" -o data/0/0-apa0.json %s' % (density, fp1, )
    print(cmd)
    os.system(cmd)

    cmd = 'wirecell-img bee-blobs -g protodunehd -s uniform -d %f --speed "1.6*mm/us" --t0 "250*us" --x0 "358*cm" -o data/0/0-apa1.json %s' % (density, fp2, )
    # cmd = 'wirecell-img bee-blobs -g protodunehd -s center --speed "1.6*mm/us" --t0 "250*us" --x0 "358*cm" -o data/0/0-apa1.json %s' % ( fp2, )
    print(cmd)
    os.system(cmd)

    cmd = 'wirecell-img bee-blobs -g protodunehd -s uniform -d %f --speed "-1.6*mm/us" --t0 "250*us" --x0 "-  358*cm" -o data/0/0-apa2.json %s' % (density, fp3, )
    print(cmd)
    os.system(cmd)

    cmd = 'wirecell-img bee-blobs -g protodunehd -s uniform -d %f --speed "1.6*mm/us" --t0 "250*us" --x0 "358*cm" -o data/0/0-apa3.json %s' % (density, fp4, )
    print(cmd)
    os.system(cmd)

    os.system('zip -r upload data')

if __name__ == "__main__":
    if (len(sys.argv)!=5):
        print("usage: python wct-img-2-bee.py 'fp1' 'fp2' 'fp3' 'fp4'")
    else:
        main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
import tarfile

file = tarfile.open(r'Downloads\jendl-f99.tar.gz')

print(file.getnames())
file.extractall(r'C:\Users\lowin\Documents\Marathon\nuclear-parser-temp\nuclear-parser-temp\My crap\JENDL fusion data')
file.close()

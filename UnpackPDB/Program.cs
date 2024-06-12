using System;
using System.IO;
using System.IO.Compression;

namespace UnpackPDB
{
    class Program
    {

        static void UnpackGzip(string fileName, string extractPath) {
            using (FileStream fInStream = new FileStream(fileName, FileMode.Open, FileAccess.Read))
            {
                using (GZipStream zipStream = new GZipStream(fInStream, CompressionMode.Decompress))
                {
                    string outName = Path.GetFileName(fileName);
                        

                    using (FileStream fOutStream = new FileStream(extractPath + outName.Remove(outName.Length - 3, 3), FileMode.Create, FileAccess.Write))
                    {
                        byte[] tempBytes = new byte[4096];
                        int i;
                        while ((i = zipStream.Read(tempBytes, 0, tempBytes.Length)) != 0)
                        {
                            fOutStream.Write(tempBytes, 0, i);
                        }
                    }
                }
            }

        }
        static void Main(string[] args)
        {
            string zipPath = Directory.GetCurrentDirectory(); 
            string gzipPath = zipPath + "\\gzips\\";      
            string extractPath = zipPath + "\\files\\";

            Directory.CreateDirectory(gzipPath);
            Directory.CreateDirectory(extractPath);

            string[] files = Directory.GetFiles(zipPath, "*.zip", SearchOption.TopDirectoryOnly);

            foreach (var file in files)
            {
                ZipFile.ExtractToDirectory(file, gzipPath);
            }

            var gzipFiles = Directory.GetFiles(gzipPath, "*.gz", SearchOption.TopDirectoryOnly);
            foreach (var file in gzipFiles)
            {
                UnpackGzip(file, extractPath);
            }
        }
    }
}

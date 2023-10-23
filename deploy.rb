#copy files into a directory structure, then zip it up.
require 'tempfile'
require 'date'

Dir.mktmpdir do |temp_dir|
  deploy_path = "./deploy/nano_recall.#{Date.today().strftime('%F')}.zip"


  deploy_i = 2
  while(File.exist?(deploy_path))
    deploy_path.gsub!(/(\.\d+)?\.zip/,".#{deploy_i}.zip")
    deploy_i += 1
  end
  deploy_path = File.expand_path(deploy_path)

  dir = "#{temp_dir}/nano_recall/"
  puts "Deploying as #{deploy_path} from #{dir}"
  FileUtils.mkdir(dir)
  FileUtils.mkdir("#{dir}/config")
  FileUtils.cp_r("lib", "#{dir}/lib")
  FileUtils.cp("config/settings.deploy.txt", "#{dir}/config/settings.txt")
  FileUtils.cp("config/reference_standards.fas", "#{dir}/config/reference_standards.fas")
  FileUtils.cp("config/HIVDB_9.0.xml", "#{dir}/config/HIVDB_9.0.xml")
  FileUtils.cp("nano_recall.rb", "#{dir}/nano_recall.rb")
  FileUtils.cp("testdata/allruns_hac/NGS059/FAP35123_HAC_NGS059_d087a0fc.fastq\ \(13\).gz", "#{dir}/example.fastq.gz")
  FileUtils.cp("README.md", "#{dir}/README.md")
  FileUtils.cp("README.md", "#{dir}/README.txt")
  FileUtils.cp("alignment_ext-0.0.0.gem", "#{dir}/alignment_ext-0.0.0.gem")
  FileUtils.cp("windows_setup.rb", "#{dir}/windows_setup.rb")

  #zip up folder and transfer it to the deploy folder.
  FileUtils.cd(temp_dir) do
    system("zip -r #{deploy_path} nano_recall")
  end
end
#folder is automatically removed at the end

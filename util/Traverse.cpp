#include <iostream>
#include "../include/UserIn.cpp"
#include <dirent.h>

#ifndef TRAVERSE_DATA
#define TRAVERSE_DATA

std::string path = "/mnt/c/";

std::vector<std::string> firstProcess(std::string directory)
{
  std::vector<std::string> f_Folders;
  std::string path_to_dir = directory;
  std::string new_path = path_to_dir + '/';
  std::vector<std::string> error1;
  error1.push_back("error... could not open directory");

  auto dir = opendir(path_to_dir.c_str());
  if (dir == NULL)
  {
    std::cout << "Could not open directory: " << directory.c_str() << std::endl;
    return error1;
  }
  auto entity = readdir(dir);
  while (entity != NULL)
  {
    if(entity->d_type == DT_DIR) 
    {
      f_Folders.push_back(std::string(entity->d_name));
    }

    entity = readdir(dir);
  }

  closedir(dir);

  return f_Folders;
}


void ProcessFile(TString file, TString path1)
{
    if (file.EndsWith(".root")) 
    {

      TFile* f;
      TTree* fTree;
      //std::cout << path1 + file << std::endl;
      f = TFile::Open((path1 + file));
      fTree = (TTree *) f->Get(global_tags::tree_tag); 
     // std::cout << "GOT TREE: " << fTree->GetName() << std::endl;
      //fTree->Print();
      if (fTree != nullptr) {
        try{
          fTree->Process(global_tags::ReaderTag); 
        } catch (...) {
          std::cout << "whoops" << std::endl;
        }
        
      }

      f->Close();
    }

}

void saveOutput(){

}

void ProcessDirectory(std::string directory, std::string current_path, std::vector<TString>* files, bool done = 0)
{
   std::string path_to_dir = current_path + directory;
   std::string new_path = path_to_dir + '/';
   std::cout << "beep" << std::endl;
   auto dir = opendir(path_to_dir.c_str());
   if (dir == NULL)
   {
      std::cout << "Could not open directory: " << directory.c_str() << std::endl;
      return;
   }

   auto entity = readdir(dir);
   while (entity != NULL)
   {
      
      if(entity->d_type == DT_DIR) 
      {
         bool do_process = true;

            if(entity->d_name[0] == '.') 
            {
               do_process = false;
            }

            if (do_process)
            {
               ProcessDirectory(std::string(entity->d_name), new_path, files, 0);
            }
        }

        if(entity->d_type == DT_REG)
        {
          //std::cout << TString(entity->d_name) << std::endl;
          TString file_name = TString(entity->d_name);
          if (file_name.EndsWith(".root")) files->push_back(TString(new_path) + TString(entity->d_name));
          //ProcessFile(std::string(entity->d_name), new_path);
        }

      entity = readdir(dir);
   }
   closedir(dir);

   if (done == 1){
      int size = files->size();
      double count = 0.0;
      double progress = size;
      for (TString file : *files){
       // std::cout << file << std::endl;
        count = count + 1.0;
        ProcessFile(file, "");
        std::cout << count / progress * 100 << "%" << std::endl;
        

      } 
   }
}


#endif
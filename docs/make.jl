push!(LOAD_PATH, "../src/");

using Documenter;
using Revise;
using MPCCLibrary;

makedocs(sitename = "MPCC Library",
        format = Documenter.HTML(
                            prettyurls = get(ENV, "CI", nothing) == "true"),
        authors = "Caroline Nakama",
        pages = ["Introduction" => "index.md",
                 "Flash Tank Problem" => "pages/flash_problem.md",
                 "Thermal Energy Storage Problem" => "pages/tes_problem.md",
                 "Bioprocess Optimization Problem" => "pages/bio_problem.md",
                 "Bilevel Optimization Problem" => "pages/biopt_problem.md",
                 "MPEC Collection" => "pages/mpec_collection.md",
                 "Methods" => "pages/mpcc_methods.md"]
        );




# flake.nix
{
  description = "A development environment for Data Science with Python and R";

  # Nix Flakes require specifying where to get packages from.
  # 'nixpkgs' is the main repository.
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  # The 'outputs' define what the flake provides.
  # We use 'flake-utils' to easily support different systems (e.g., x86_64, aarch64).
  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        # Import the package set for the specific system.
        pkgs = import nixpkgs {
          inherit system;
          config.allowUnfree = true; # Allow packages with unfree licenses if needed.
        };

        # Define Python packages, translating from your conda list.
        python-packages = ps: with ps; [
          # Python interpreter and base packages
          python311
          # Core data science
          pandas
          numpy
          scipy
          scikit-learn
          # From your conda list
          scanpy
          leidenalg
          # Radian shell for R
          radian
        ];

        # Create a Python environment with the specified packages.
        python-env = pkgs.python311.withPackages python-packages;

      in
      {
        # This is the main development shell, equivalent to your devcontainer.
        devShells.default = pkgs.mkShell {
          name = "datascience-shell";

          # List of packages available in the shell, from your apt-get list.
          packages = with pkgs; [
            # Base development tools
            coreutils # for basic commands
            bashInteractive
            gcc
            gfortran
            git
            gh # GitHub CLI
            htop
            vim
            parallel
            wget
            curl
            gnumake

            # R environment
            R
            r-cran-docopt # From your apt-get list
            
            # Python Environment defined above
            python-env

            # System libraries required for R/Python packages, from your apt-get list.
            # Nix automatically manages dependencies, so you often need fewer explicit libraries.
            openssl
            libcurl
            libxml2
            zlib
            hdf5 # For hdf5r and h5py
            gmp # libgmp-dev
            openblas # libopenblas-dev
            tbb # libtbb-dev
            fontconfig # libfontconfig1
            cairo
            pango

            # Quarto CLI
            quarto

            # Node.js and npm
            nodejs
            pkgs.nil  # The language server
            pkgs.nixpkgs-fmt # The formatter (often preferred over nixfmt for its wider adoption)
          ];

          # This is the equivalent of your postCreateCommand.sh and .Rprofile setup.
          # It runs every time you enter the environment.
          shellHook = ''
            echo "✅ Data Science environment activated."

            # Set environment variables, similar to Dockerfile ENV.
            export LANG="en_US.UTF-8"
            export LC_ALL="en_US.UTF-8"
            export TZ="UTC"
            # The PATH is automatically configured by Nix to include all packages.

            # Set git config, as in your postCreateCommand.sh
            git config --global user.name "h4rvey-g"
            git config --global user.email "babaolanqiu@gmail.com"

            # Create a default .Rprofile if it doesn't exist
            # This is a better practice than appending, as it's idempotent.
            if [ ! -f ~/.Rprofile ]; then
              echo 'options(defaultPackages=c(getOption("defaultPackages"), "tidyverse", "targets", "skimr", "gittargets"))' > ~/.Rprofile
              echo 'if (interactive() && Sys.getenv("RSTUDIO") == "") {
                source(file.path(Sys.getenv(if (.Platform$OS.type == "windows") "USERPROFILE" else "HOME"), ".vscode-R", "init.R"))
              }' >> ~/.Rprofile
            fi
            
            echo "💡 Type 'R' to start R, or 'radian' for a better R console."
            echo "💡 Type 'python' to start Python."
          '';
        };
      });
}

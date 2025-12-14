
"""
HiFAST Global Configuration Module
----------------------------------

This module handles global configuration management for the HiFAST package.
It supports layered configuration with the following priority (highest first):

1. **Environment Variables**: `HIFAST_<SECTION>_<OPTION>` (e.g., `HIFAST_TCAL_SOURCE`)
2. **User Config File**: `~/.config/hifast/config.ini` (standard INI format)
3. **Defaults**: Hardcoded defaults within this module.

Usage for Developers
--------------------
To add configuration for a new module (e.g., 'fitting'):

1. **Register the Option**:
   Add your new option to `_KNOWN_OPTIONS` dictionary below. This serves as
   direct documentation and enables the `list` CLI command.
   
   _KNOWN_OPTIONS = {
       ...
       'fitting': {
           'max_iter': {'default': 50, 'help': "Maximum iterations for fitter"},
           'algorithm': {'default': 'gaussian', 'help': "Fitting algorithm to use"}
       }
   }

2. **Access in Code**:
   Import `conf` and use typed getters:
   
   from hifast.utils.config import conf
   max_iter = conf.get_int('fitting', 'max_iter', 50)
   algo = conf.get('fitting', 'algorithm', 'gaussian')

"""

import os
import sys
import configparser
from pathlib import Path

# Constants for Defaults
_TCAL_GITHUB_MANIFEST = "https://raw.githubusercontent.com/jyingjie/hifast-tcal-data/main/manifest.json"
_TCAL_GITHUB_BASE = "https://github.com/jyingjie/hifast-tcal-data/releases/download"

_TCAL_R2_MANIFEST = "https://pub-b2347aee569d45789e90c47ad3177e0e.r2.dev/tcal-data/manifest.json"
_TCAL_R2_BASE = "https://pub-b2347aee569d45789e90c47ad3177e0e.r2.dev/tcal-data"

# Known Options Registry
# ----------------------
# DOCUMENTATION SOURCE OF TRUTH.
# Add new module parameters here to expose them to users via `hifast.utils.config list`.
_KNOWN_OPTIONS = {
    'general': {
        'offline': {
            'default': False,
            'help': "Enable offline mode. Disable all network requests. (Env: HIFAST_GENERAL_OFFLINE)"
        }
    },
    'tcal': {
        'source': {
            'default': 'r2', 
            'help': "Data source used for Tcal files. Options: 'r2' (default, Cloudflare) or 'github'."
        },
        'tcal_dir': {
            'default': '~/Tcal/', 
            'help': "Local directory where Tcal data is stored/cached."
        },
        'manifest_url': {
            'default': None, 
            'help': "Remote Manifest URL. Required if source='custom'. Can also be used to override 'r2'/'github' defaults."
        },
        'base_url': {
            'default': None, 
            'help': "Remote Base URL. Required if source='custom'. Can also be used to override 'r2'/'github' defaults."
        },
    }
}

_TCAL_SOURCES = {
    'github': {'manifest': _TCAL_GITHUB_MANIFEST, 'base': _TCAL_GITHUB_BASE},
    'r2': {'manifest': _TCAL_R2_MANIFEST, 'base': _TCAL_R2_BASE}
}

# XDG Config Home or ~/.config
CONFIG_DIR = Path(os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config")) / "hifast"
CONFIG_FILE = CONFIG_DIR / "config.ini"

class ConfigManager:
    def __init__(self):
        self._config = configparser.ConfigParser()
        self._load()

    def _load(self):
        """Loads configuration from file."""
        if CONFIG_FILE.exists():
            try:
                self._config.read(CONFIG_FILE)
            except Exception as e:
                print(f"Warning: Failed to load config file {CONFIG_FILE}: {e}", file=sys.stderr)
        
        # Ensure default sections exist
        if 'tcal' not in self._config:
            self._config['tcal'] = {}
        if 'general' not in self._config:
            self._config['general'] = {}

    def save(self):
        """Saves current configuration to file."""
        try:
            CONFIG_DIR.mkdir(parents=True, exist_ok=True)
            with open(CONFIG_FILE, 'w') as f:
                self._config.write(f)
        except Exception as e:
            print(f"Error saving config file: {e}", file=sys.stderr)

    def get(self, section, option, default=None):
        """
        Get value with priority:
        1. Environment Variable (HIFAST_SECTION_OPTION)
        2. Config File
        3. Default
        """

        env_var = f"HIFAST_{section.upper()}_{option.upper()}"
        if env_var in os.environ:
            return os.environ[env_var]
        
        return self._config.get(section, option, fallback=default)

    def get_int(self, section, option, default=None):
        """Get integer value with environment override."""
        val = self.get(section, option, default)
        try:
            return int(val)
        except (ValueError, TypeError):
            return default

    def get_float(self, section, option, default=None):
        """Get float value with environment override."""
        val = self.get(section, option, default)
        try:
            return float(val)
        except (ValueError, TypeError):
            return default

    def get_boolean(self, section, option, default=None):
        """Get boolean value with environment override."""
        val = self.get(section, option, default)
        if val is None:
            return default
        if isinstance(val, bool):
            return val
        s = str(val).lower()
        if s in ('true', 'yes', 'on', '1'):
            return True
        if s in ('false', 'no', 'off', '0'):
            return False
        return default

    def set(self, section, option, value, save=False):
        """Sets a value in memory, optionally saving to file."""
        if section not in self._config:
            self._config[section] = {}
        self._config[section][option] = str(value)
        if save:
            self.save()

    def get_tcal_urls(self):
        """
        Resolves Tcal URLs based on configuration.
        Returns: (manifest_url, base_url)
        """
        # 1. Determine base defaults from source
        source = self.get('tcal', 'source', 'r2')
        
        preset = _TCAL_SOURCES.get(source, {})
        default_manifest = preset.get('manifest')
        default_base = preset.get('base')
        
        if not default_manifest and source != 'custom':
             print(f"Warning: Unknown [tcal] source='{source}'. Using GitHub default.", file=sys.stderr)
             default_manifest = _TCAL_GITHUB_MANIFEST
             default_base = _TCAL_GITHUB_BASE

        # 2. Check for explicit overrides (Env or Config)
        # These take precedence over source-derived defaults
        manifest_url = self.get('tcal', 'manifest_url')
        base_url = self.get('tcal', 'base_url')

        return (manifest_url or default_manifest), (base_url or default_base)

# Global instance
# Global instance
conf = ConfigManager()

def get_parser():
    import argparse
    import textwrap
    
    description = textwrap.dedent("""
    Manage HiFAST Global Configuration
    ==================================
    
    Set and retrieve configuration values for HiFAST modules (e.g., Tcal).
    
    Priority Order:
      1. Environment Variables:
         - Standard: HIFAST_SECTION_OPTION (e.g., HIFAST_TCAL_SOURCE)
      2. Config File (~/.config/hifast/config.ini)
      3. Global Defaults
    """)
    
    epilog = textwrap.dedent("""
    Examples:
      # List all options
      python -m hifast.utils.config list
      
      # Set Tcal source to R2 (default)
      python -m hifast.utils.config set tcal source r2
      
      # Get current data directory
      python -m hifast.utils.config get tcal data_dir
    """)

    parser = argparse.ArgumentParser(
        description=description, 
        epilog=epilog, 
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparsers = parser.add_subparsers(dest='command')
    
    # set command
    set_parser = subparsers.add_parser('set', help='Set a configuration value')
    set_parser.add_argument('section', help='Config section (e.g., tcal)')
    set_parser.add_argument('option', help='Config option (e.g., source)')
    set_parser.add_argument('value', help='Value to set')
    
    # get command
    get_parser_cmd = subparsers.add_parser('get', help='Get a configuration value')
    get_parser_cmd.add_argument('section', help='Config section')
    get_parser_cmd.add_argument('option', help='Config option')

    # list command
    subparsers.add_parser('list', help='List all available options and current values')
    
    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    
    if args.command == 'set':
        # Validation checks
        if args.section in _KNOWN_OPTIONS:
            valid_opts = _KNOWN_OPTIONS[args.section].keys()
            if args.option not in valid_opts:
                import sys
                import difflib
                print(f"Warning: '{args.option}' is not a known option for section '{args.section}'.", file=sys.stderr)
                matches = difflib.get_close_matches(args.option, valid_opts)
                if matches:
                    print(f"Did you mean: {', '.join(matches)}?", file=sys.stderr)
                print(f"Valid options: {', '.join(valid_opts)}", file=sys.stderr)
        else:
             import sys
             print(f"Warning: '{args.section}' is not a known section.", file=sys.stderr)
             print(f"Known sections: {', '.join(_KNOWN_OPTIONS.keys())}", file=sys.stderr)

        conf.set(args.section, args.option, args.value, save=True)
        print(f"Set [{args.section}] {args.option} = {args.value}")
    elif args.command == 'get':
        # Lookup default from KNOWN_OPTIONS if available
        default = None
        if args.section in _KNOWN_OPTIONS and args.option in _KNOWN_OPTIONS[args.section]:
            default = _KNOWN_OPTIONS[args.section][args.option].get('default')
            
        val = conf.get(args.section, args.option, default)
        print(val)
    elif args.command == 'list':
        import textwrap
        col_sec = 10
        col_opt = 15
        col_env = 30
        col_val = 30
        col_desc = 35 # Allow simple wrap
        
        # Header
        print(f"{'SECTION':<{col_sec}} {'OPTION':<{col_opt}} {'ENV VAR':<{col_env}} {'VALUE':<{col_val}} {'DESCRIPTION'}")
        print("-" * 130)
        
        # 1. Print Known Options
        for section, options in _KNOWN_OPTIONS.items():
            for opt, meta in options.items():
                val = str(conf.get(section, opt, meta.get('default')))
                desc = meta.get('help', '')
                env_var = f"HIFAST_{section.upper()}_{opt.upper()}"
                
                # Truncate value if too long (rare config)
                if len(val) > col_val - 2:
                    val = val[:col_val-5] + "..."
                
                # Wrap description
                wrapped_desc = textwrap.wrap(desc, width=50) # Use loose width for desc
                
                # First line
                d_first = wrapped_desc[0] if wrapped_desc else ""
                print(f"{section:<{col_sec}} {opt:<{col_opt}} {env_var:<{col_env}} {val:<{col_val}} {d_first}")
                
                # Subsequent lines
                for line in wrapped_desc[1:]:
                     print(f"{'':<{col_sec}} {'':<{col_opt}} {'':<{col_env}} {'':<{col_val}} {line}")
    else:
        parser.print_help()

if __name__ == "__main__":
    main()

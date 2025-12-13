
import os
import sys
import unittest
import shutil
import tempfile
from pathlib import Path
from unittest.mock import patch

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

class TestConfigSystem(unittest.TestCase):
    def setUp(self):
        # Use a temporary directory for config file
        self.test_dir = tempfile.mkdtemp()
        self.config_dir = Path(self.test_dir) / "hifast"
        
        # Patch XDG_CONFIG_HOME in os.environ so reload picks it up
        self.env_patcher = patch.dict(os.environ, {"XDG_CONFIG_HOME": self.test_dir})
        self.env_patcher.start()
        
        # Reload config to re-calculate CONFIG_DIR/CONFIG_FILE based on new env
        from hifast.utils import config
        import importlib
        importlib.reload(config)
        self.config_module = config
        self.conf = config.conf

    def tearDown(self):
        self.env_patcher.stop()
        shutil.rmtree(self.test_dir)

    def test_default_urls_is_r2(self):
        """Test default is now R2"""
        m, b = self.conf.get_tcal_urls()
        self.assertIn("pub-b2347", m)

    def test_tcal_dir(self):
        """Test tcal tcal_dir configuration"""
        # Default
        self.assertEqual(self.conf.get('tcal', 'tcal_dir', '~/Tcal/'), '~/Tcal/')
        
        # Override
        self.conf.set('tcal', 'tcal_dir', '/custom/path')
        self.assertEqual(self.conf.get('tcal', 'tcal_dir'), '/custom/path')

    def test_typed_getters(self):
        """Test get_int, get_float, get_boolean"""
        self.conf.set('test', 'myint', '42')
        self.conf.set('test', 'myfloat', '3.14')
        self.conf.set('test', 'mybool', 'true')
        
        self.assertEqual(self.conf.get_int('test', 'myint'), 42)
        self.assertEqual(self.conf.get_float('test', 'myfloat'), 3.14)
        self.assertTrue(self.conf.get_boolean('test', 'mybool'))
        
        # Test defaults
        self.assertEqual(self.conf.get_int('test', 'missing', 10), 10)
        self.assertFalse(self.conf.get_boolean('test', 'missing', False))

if __name__ == "__main__":
    unittest.main()

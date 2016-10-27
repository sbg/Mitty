import pkg_resources

# These test files need to be available to the rest of the test suite
# They are included in MANIFEST.in
example_data_dir = pkg_resources.resource_filename(__name__, 'data')
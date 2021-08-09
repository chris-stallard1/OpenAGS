from setuptools import setup
setup(
    name="openags",
    version='0.1.7',
    packages=['openags'],
    install_requires=['scipy','numpy','sigfig','xylib-py-wheels'],
    classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Physics',
          'Topic :: Software Development :: Libraries :: Python Modules',
          ],
      author='Christopher Stallard',
      author_email='christopher.stallard1@gmail.com',
      license='MIT',
      url='https://github.com/chris-stallard1/OpenAGS'
      )

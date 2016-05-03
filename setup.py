from setuptools import setup
import pip.req

install_reqs = pip.req.parse_requirements('requirements.txt', session=False)
reqs = [str(ir.req) for ir in install_reqs]

setup(name='PDS',
      version='1.0',
      description='Poisson disk sampling approach for brain parcellations.',
      url='https://github.com/mdschirmer/poisson_disk_sampling.git',
      author='Markus D. Schirmer',
      author_email='software@markus-schirmer.com',
      license='MIT',
      packages=['PDS'],
      install_requires=reqs,
      zip_safe=False)
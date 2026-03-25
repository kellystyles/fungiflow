from setuptools import setup, find_packages

setup(
    name='fungiflow',
    version='0.1',
    packages=find_packages(),
    install_requires=[],  # Add your dependencies here
    author='Kelly Styles',
    author_email='your_email@example.com',  # Update with your email
    description='A brief description of your package.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/kellystyles/fungiflow',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
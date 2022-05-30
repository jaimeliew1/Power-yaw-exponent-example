# Power-yaw-exponent-example

A recreation of the analysis in the journal paper titled *Analytical model for the power–yaw sensitivity of wind turbines operating in full wake* ([doi.org/10.5194/wes-5-427-2020](http://doi.org/10.5194/wes-5-427-2020)).

## Installation
1) Clone this repository and use the package manager [pip](https://pip.pypa.io/en/stable/) to install the dependencies (which are numpy, scipy, and matplotlib).

    ```bash
    git clone https://github.com/jaimeliew1/Power-yaw-exponent-example.git
    pip install foobar
    ```

2) Download the data set *LES of wake flow behind 2.3MW wind turbine* ([https://doi.org/10.11583/DTU.12005421.v1](https://doi.org/10.11583/DTU.12005421.v1)) and place the `.*fluc` files in the `data` folder of this repository.



## Usage
Run the `process_dataset.py`script in Python:
```bash
python src/process_dataset.py
```
This should produce a table of alpha values for a given downstream distance in the `out`folder.
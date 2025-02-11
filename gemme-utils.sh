# Set GU_PATH to this directory
export GU_PATH=$(readlink -f "$(dirname "${BASH_SOURCE[0]}")")


# Function to execute Python scripts from GU_PATH
# Usage: gemme-utils <script.py> [args...]
gemme-utils() {

    if [[ -z "$1" ]]; then
        echo "Usage: gemme-utils <script.py> [args...]"
        return 1
    elif [[ ! -f "$GU_PATH/$1" ]]; then
        echo "Error: Script '$1' not found"
        return 1
    fi

    python3 "$GU_PATH/$1" "${@:2}"
}


# Bash Completion Function for gemme-utils
_gemme_utils_complete() {

    _complete_scripts() {
        local cur=${COMP_WORDS[COMP_CWORD]}
        local scripts

        mapfile -t scripts < <(find "$GU_PATH" -maxdepth 1 -type f -name "*.py" -exec basename {} \;)
        
        COMPREPLY=( $(compgen -W "${scripts[*]}" -- "$cur") )
    }

    _complete_default() {
        local script=${COMP_WORDS[1]}  
        local cur=${COMP_WORDS[COMP_CWORD]} 

        if [ -f "$GU_PATH/$script" ]; then
            COMPREPLY=( $(compgen -f -- $cur) )  # Default suggestion
        fi
    }

    if [ ${COMP_CWORD} -eq 1 ]; then
        _complete_scripts
    else
        _complete_default
    fi
}

complete -F _gemme_utils_complete gemme-utils

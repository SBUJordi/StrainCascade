#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# straincascade_theme.sh
# CLI theming, colors, logo, and display formatting functions

## Color definitions (portable ANSI codes)
# Check if terminal supports colors
if [[ -t 1 ]] && [[ -n "$TERM" ]] && command -v tput &>/dev/null && [[ $(tput colors 2>/dev/null) -ge 8 ]]; then
    readonly COLOR_ENABLED=1
else
    readonly COLOR_ENABLED=0
fi

# Define colors
if [[ "$COLOR_ENABLED" -eq 1 ]]; then
    readonly C_RESET='\e[0m'
    readonly C_BOLD='\e[1m'
    readonly C_DIM='\e[2m'
    # Primary: Orange accent
    readonly C_ORANGE='\e[38;2;250;128;0m'
    # Secondary: Light blue accent  
    readonly C_BLUE='\e[38;2;135;206;250m'
    # Text colors
    readonly C_WHITE='\e[97m'
    readonly C_GRAY='\e[90m'
    readonly C_GREEN='\e[38;2;80;200;120m'
    readonly C_YELLOW='\e[38;2;255;200;87m'
    readonly C_RED='\e[38;2;255;99;71m'
else
    readonly C_RESET=''
    readonly C_BOLD=''
    readonly C_DIM=''
    readonly C_ORANGE=''
    readonly C_BLUE=''
    readonly C_WHITE=''
    readonly C_GRAY=''
    readonly C_GREEN=''
    readonly C_YELLOW=''
    readonly C_RED=''
fi

## Logo and branding functions
print_logo() {
    local max_height heights colors widths aligned_bottom
    local line row idx h color_type w print_bar local_k start_line
    local intensity r g b full_r full_g full_b
    
    # Function to print a row of bars with gradient
    _print_bars() {
        local max_height=$1
        local aligned_bottom=${!#}
        local -a heights colors widths
        for i in {2..14}; do heights+=("${!i}"); done
        for i in {15..27}; do colors+=("${!i}"); done  
        for i in {28..40}; do widths+=("${!i}"); done
        
        for line in $(seq 1 "$max_height"); do
            row=""
            for idx in "${!heights[@]}"; do
                [[ $idx -gt 0 ]] && row="$row "
                h=${heights[$idx]}
                color_type=${colors[$idx]}
                w=${widths[$idx]}
                print_bar=0
                local_k=0
                
                if [[ "$aligned_bottom" -eq 1 ]]; then
                    start_line=$((max_height - h + 1))
                    if [[ "$line" -ge "$start_line" ]]; then
                        print_bar=1
                        local_k=$((line - start_line + 1))
                    fi
                else
                    if [[ "$line" -le "$h" ]]; then
                        print_bar=1
                        local_k=$line
                    fi
                fi
                
                if [[ "$print_bar" -eq 1 ]]; then
                    if [[ "$color_type" = "O" ]]; then
                        full_r=250; full_g=128; full_b=0
                    else
                        full_r=135; full_g=206; full_b=250
                    fi
                    
                    if [[ "$aligned_bottom" -eq 1 ]]; then
                        intensity=$(( local_k * 100 / h ))
                    else
                        intensity=$(( (h - local_k + 1) * 100 / h ))
                    fi
                    
                    r=$(( intensity * full_r / 100 ))
                    g=$(( intensity * full_g / 100 ))
                    b=$(( intensity * full_b / 100 ))
                    
                    if [[ "$COLOR_ENABLED" -eq 1 ]]; then
                        row="${row}\e[38;2;${r};${g};${b}m$(printf '█%.0s' $(seq 1 $w))\e[0m"
                    else
                        row="${row}$(printf '█%.0s' $(seq 1 $w))"
                    fi
                else
                    row="${row}$(printf ' %.0s' $(seq 1 $w))"
                fi
            done
            echo -e "$row"
        done
    }
    
    local widths_arr=(2 2 2 2 2 2 1 1 1 1 1 1 1)
    local upper_heights=(8 5 7 2 5 3 2 4 3 6 2 4 2)
    local upper_colors=(O O O B O B O B B B O B B)
    local lower_heights=(2 3 2 5 3 4 5 4 6 3 7 5 8)
    local lower_colors=(O O O B O B O B B B O B B)
    
    echo
    _print_bars 8 "${upper_heights[@]}" "${upper_colors[@]}" "${widths_arr[@]}" 1
    echo -e "         ${C_BOLD}${C_GRAY}StrainCascade${C_RESET}"
    _print_bars 8 "${lower_heights[@]}" "${lower_colors[@]}" "${widths_arr[@]}" 0
    echo -e "${C_DIM}${C_GRAY}by Sebastian Bruno Ulrich JORDI${C_RESET}"
    echo
}

print_header() {
    print_logo
}

## Section and formatting functions
print_section() {
    local title="$1"
    echo -e "\n${C_BOLD}${C_BLUE}${title}${C_RESET}"
}

print_option() {
    local opt="$1"
    local desc="$2"
    printf "  ${C_BLUE}%-28s${C_RESET} %s\n" "$opt" "$desc"
}

print_sub() {
    local text="$1"
    echo -e "  ${C_DIM}                             ${text}${C_RESET}"
}

print_example() {
    local desc="$1"
    local cmd="$2"
    echo -e "  ${C_GRAY}# ${desc}${C_RESET}"
    echo -e "  ${C_GREEN}${cmd}${C_RESET}"
    echo
}

print_module() {
    local id="$1"
    local name="$2"
    local desc="$3"
    printf "  ${C_WHITE}%-4s${C_RESET} %-36s ${C_DIM}%s${C_RESET}\n" "$id" "$name" "$desc"
}

## Message helpers for consistent CLI output
print_error() {
    echo -e "${C_RED}Error:${C_RESET} $1" >&2
}

print_warning() {
    echo -e "${C_YELLOW}Warning:${C_RESET} $1" >&2
}

print_success() {
    echo -e "${C_GREEN}✓${C_RESET} $1"
}

print_info() {
    echo -e "${C_BLUE}ℹ${C_RESET} $1"
}

print_status() {
    local label="$1"
    local value="$2"
    printf "  ${C_BLUE}%-14s${C_RESET} %s\n" "${label}:" "$value"
}

print_divider() {
    echo -e "${C_DIM}─────────────────────────────────────────────────────────────────${C_RESET}"
}

import { Box, Select, MenuItem, FormControl } from '@mui/material'
import { useEffect, useState } from 'react'

const currencies = [
    {
        value: 'USD',
        label: '$',
    },
    {
        value: 'EUR',
        label: '€',
    },
    {
        value: 'BTC',
        label: '฿',
    },
    {
        value: 'JPY',
        label: '¥',
    },
];

export function InlineSelector(){

    const [currency, setCurrency] = useState('EUR');

    const handleChange = (event) => {
        setCurrency(event.target.value);
    };

    return (
        <FormControl variant="standard" disabled>
            <Select
            value={currency}
            onChange={handleChange}
            label="Age"
            >
            {currencies.map((option) => (
                <MenuItem key={option.value} value={option.value}>
                {option.label}
                </MenuItem>
            ))}
            </Select>
        </FormControl>
    )
}

export function FilterItem({variant, content}){

    useEffect(()=>{

    }, [])

    return (
        <Box>
            
        </Box>
    )
}

export function Filter({children, width="100%", height="100%"}){
    
    return (
        <Box
            sx={{
                width, height, 
                display: "flex",
                flexDirection: "row",
                justifyContent: "space-between",
                padding: "3%",
                backgroundColor: "white",
                borderRadius: "20px",
                boxShadow: "0px 4px 6px rgba(0, 0, 0, 0.15), 0px 0px 1px rgba(0, 0, 0, 0.4)",
                zIndex: "100"
            }}
        >
            {children}
        </Box>
    )
}
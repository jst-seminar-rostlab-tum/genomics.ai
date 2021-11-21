import React from 'react';
import { createTheme, ThemeProvider, styled } from '@mui/material/styles';
import Switch from '@mui/material/Switch';

const SwitchButton = styled((props) => (
  <Switch focusVisibleClassName=".Mui-focusVisible" disableRipple {...props} />
))(({ theme }) => ({
  width: 42,
  height: 26,
  padding: 0,
  '& .MuiSwitch-switchBase': {
    padding: 0,
    margin: 2,
    transitionDuration: '300ms',
    '&.Mui-checked': {
      transform: 'translateX(16px)',
      color: '#fff',
      '& + .MuiSwitch-track': {
        backgroundColor: theme.palette.mode === 'dark' ? '#1888ff' : '#1888ff',
        opacity: 1,
        border: 0,
      },
      '&.Mui-disabled + .MuiSwitch-track': {
        opacity: 0.5,
      },
    },
    '&.Mui-focusVisible .MuiSwitch-thumb': {
      color: '#1888ff',
      border: '6px solid #fff',
    },
    '&.Mui-disabled .MuiSwitch-thumb': {
      color:
                theme.palette.mode === 'light'
                  ? theme.palette.grey[100]
                  : theme.palette.grey[600],
    },
    '&.Mui-disabled + .MuiSwitch-track': {
      opacity: theme.palette.mode === 'light' ? 0.7 : 0.3,
    },
  },
  '& .MuiSwitch-thumb': {
    boxSizing: 'border-box',
    width: 22,
    height: 22,
  },
  '& .MuiSwitch-track': {
    borderRadius: 26 / 2,
    backgroundColor: theme.palette.mode === 'light' ? 'darkgray' : 'darkgray',
    opacity: 1,
    transition: theme.transitions.create(['background-color'], {
      duration: 500,
    }),
  },
}));
export default SwitchButton;

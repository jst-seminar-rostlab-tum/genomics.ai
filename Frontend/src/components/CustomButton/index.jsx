import * as React from 'react';
import { ButtonBase } from '@mui/material';
import { colors } from 'shared/theme/colors';
import styles from './custombutton.module.css';

function getStyles(type ) {
  switch (type) {
    case 'secondary':
      return ({
        backgroundColor: colors.secondary1[400],
        borderRadius: "100px",
        padding: "0.5em 1em 0.5em 1em",
        fontSize: "1em",
        color: 'white',
        ':hover': { backgroundColor: colors.secondary1[500], transition: '0.4s', color: 'white', boxShadow: "0px 4px 20px rgba(160, 230, 158, 0.2)" },
        ':active': { backgroundColor: `${colors.secondary1[600]} !important`, transition: '0.4s', color: 'white' },
        ':focus': { backgroundColor: colors.secondary1[300], transition: '0.4s', color: 'white' },
        ':disabled': { backgroundColor: '#DAF3DB', transition: '0.4s', color: colors.secondary1[600] },
      });
    case 'tertiary':
      return ({
        color: colors.neutral[900],
        borderRadius: "100px",
        padding: "0.5em 1em 0.5em 1em",
        fontSize: "1em",
        background: "none",
        ':hover': { transition: '0.4s', color: colors.primary[500] },
        ':focus': { transition: '0.4s', color: colors.primary[500], textDecoration: 'underline' },
        ':disabled': { transition: '0.4s', color: colors.neutral[400] },
      });
    default:
      return ({
        backgroundColor: colors.primary[400],
        borderRadius: "100px",
        padding: "0.5em 1em 0.5em 1em",
        fontSize: "1em",
        color: 'white',
        ':hover': { backgroundColor: colors.primary[500], transition: '0.4s', color: 'white', boxShadow: "0px 4px 20px rgba(24, 64, 96, 0.2)" },
        ':active': { backgroundColor: `${colors.primary[600]} !important`, transition: '0.4s', color: 'white' },
        ':focus': { backgroundColor: colors.primary[300], transition: '0.4s', color: 'white' },
        ':disabled': { backgroundColor: '#EBEFFF', transition: '0.4s', color: colors.primary[600] },
      });
  }
}


/**
 * 
 * Custom Button 
 * Accepts all MUI props
 * @param size ['small', 'medium', 'large']
 * @param type ['primary', 'secondary', 'tertiary']
 * https://www.figma.com/file/HcTwUyNxjZJksJ7yfncM9y/JST-22-Design?node-id=42%3A2
 */
const CustomButton = (props) => { 

  console.log(props.size)
  return (
  <ButtonBase
    disableRipple
    {...props}
    sx={getStyles(props.type, props.size)}
    className={styles.button}
    variant={props.type === 'tertiary' ? 'text' : 'contained'}
  >
    {props.children}
  </ButtonBase>
  )
};

export default CustomButton;

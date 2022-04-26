import * as React from 'react';
import Button from '@mui/material/Button';
import { colors } from 'shared/theme/colors';
import styles from './custombutton.module.css';

/**
 * Customized button
 * @param {string} text - The text that is displayed inside the button
 * @param {function} onClick - The click handler function
 * @param {boolean} disabled
 * @param {string} size - 'small', 'medium' or 'large', default is 'medium'
 * @param {string} type - 'primary', 'secondary' or 'tertiary', default is 'primary'
 * https://www.figma.com/file/HcTwUyNxjZJksJ7yfncM9y/JST-22-Design?node-id=42%3A2
 */

/* Example usage
Import CustomButton to any component, see below examples:
<CustomButton text="Tertiary" onClick={handleButtonClick} type="tertiary"/>
<CustomButton text="Secondary" onClick={handleButtonClick} type="secondary" />
<CustomButton text="DisabledLargePrimary" onClick={handleButtonClick} size="large" disabled />
*/

function getStyles(type) {
  switch (type) {
    case 'secondary':
      return ({
        backgroundColor: colors.secondary1[400],
        color: 'white',
        '&:hover': { backgroundColor: colors.secondary1[500], transition: '0.4s', color: 'white' },
        '&:active': { backgroundColor: `${colors.secondary1[600]} !important`, transition: '0.4s', color: 'white' },
        '&:focus': { backgroundColor: colors.secondary1[300], transition: '0.4s', color: 'white' },
        '&:disabled': { backgroundColor: '#DAF3DB', transition: '0.4s', color: colors.secondary1[600] },
      });
    case 'tertiary':
      return ({
        color: colors.neutral[900],
        '&:hover': { transition: '0.4s', color: colors.primary[500] },
        '&:focus': { transition: '0.4s', color: colors.primary[500], textDecoration: 'underline' },
        '&:disabled': { transition: '0.4s', color: colors.neutral[400] },
      });
    default:
      return ({
        backgroundColor: colors.primary[400],
        color: 'white',
        '&:hover': { backgroundColor: colors.primary[500], transition: '0.4s', color: 'white' },
        '&:active': { backgroundColor: `${colors.primary[600]} !important`, transition: '0.4s', color: 'white' },
        '&:focus': { backgroundColor: colors.primary[300], transition: '0.4s', color: 'white' },
        '&:disabled': { backgroundColor: '#EBEFFF', transition: '0.4s', color: colors.primary[600] },
      });
  }
}

const CustomButton = ({
  text, onClick, disabled = false, size = 'medium', type = 'primary',
}) => (
  <Button
    disableRipple
    size={size}
    sx={getStyles(type)}
    onClick={onClick}
    disabled={disabled}
    className={styles.button}
    variant={type === 'tertiary' ? 'text' : 'contained'}
  >
    {text}
  </Button>
);

export default CustomButton;

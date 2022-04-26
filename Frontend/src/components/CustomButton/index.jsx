import * as React from 'react';
import Button from '@mui/material/Button';
import { colors } from 'shared/theme/colors';

/* Example usage
Import CustomButton to any component, see below examples:
<CustomButton text="Tertiary" onClick={handleButtonClick} type="tertiary"/>
<CustomButton text="Secondary" onClick={handleButtonClick} type="secondary" />
<CustomButton text="DisabledLargePrimary" onClick={handleButtonClick} size="large" disabled=true />
 */

function getStyles(type) {
  switch (type) {
    case 'secondary':
      return ({
        backgroundColor: colors.secondary1[400],
        color: 'white',
        margin: 2,
        textTransform: 'none',
        '&:hover': { backgroundColor: colors.secondary1[500], transition: '0.4s', color: 'white' },
        '&:active': { backgroundColor: `${colors.secondary1[600]} !important`, transition: '0.4s', color: 'white' },
        '&:focus': { backgroundColor: colors.secondary1[300], transition: '0.4s', color: 'white' },
        '&:disabled': { backgroundColor: '#DAF3DB', transition: '0.4s', color: colors.secondary1[600] },
      });
    case 'tertiary':
      return ({
        color: colors.neutral[900],
        margin: 2,
        textTransform: 'none',
        '&:hover': { transition: '0.4s', color: colors.primary[500] },
        '&:focus': { transition: '0.4s', color: colors.primary[500], textDecoration: 'underline' },
        '&:disabled': { transition: '0.4s', color: colors.neutral[400] },
      });
    default:
      return ({
        backgroundColor: colors.primary[400],
        color: 'white',
        margin: 2,
        textTransform: 'none',
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
    variant={type === 'tertiary' ? 'text' : 'contained'}
  >
    {text}
  </Button>
);

export default CustomButton;

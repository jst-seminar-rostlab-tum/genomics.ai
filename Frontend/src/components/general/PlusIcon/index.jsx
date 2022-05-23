import {
  Icon, IconButton, createTheme, ThemeProvider,
} from '@mui/material';
import React from 'react';
import { colors } from 'shared/theme/colors';
import AddIcon from '@mui/icons-material/Add';

// const themeIcon = createTheme({
//   palette: {
//     primary: colors.plus,
//   },
// });
// use primary as backgroundColor in line 21 not working. So I just hard coded the color
// in order to let the plus thinner,we can not use the mui plusicon.

export default function PlusIcon({ onClick }) {
  return (
  // <ThemeProvider theme={themeIcon}>
    <IconButton
      sx={{ p: 0 }}
      style={{
        backgroundColor: '#5676E4', minHeight: 38, minWidth: 38,
      }}
      aria-label="plus icon"
      onClick={onClick}
    >
      <AddIcon sx={{ fontSize: 30, color: 'white' }} />
    </IconButton>
  // </ThemeProvider>

  );
}

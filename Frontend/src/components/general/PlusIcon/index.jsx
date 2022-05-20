import {
  Icon, IconButton, createTheme, ThemeProvider
} from '@mui/material';
import React from 'react';
import { colors } from 'shared/theme/colors';

const themeIcon = createTheme({
  palette: {
    primary: colors.plus,
  },
});

export default function PlusIcon({ onClick }) {
  return (
    <ThemeProvider theme={themeIcon}>
      <IconButton sx={{ p: 0 }} aria-label="plus icon" onClick={onClick}>
        <Icon
          className="plusicon"
          color="primary"
          sx={{ fontSize: 45 }}
        >
          add_circle

        </Icon>
      </IconButton>
    </ThemeProvider>

  );
}

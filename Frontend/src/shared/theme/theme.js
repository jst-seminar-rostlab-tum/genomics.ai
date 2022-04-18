import { createTheme } from "@mui/material/styles" 
import { colors } from "./colors"

export const theme = createTheme({
  palette: {
    primary: {
      dark: colors.primary[900],
      main: colors.primary[700],
      light: colors.primary[400],
      contrastText: colors.neutral[100]
    },
    secondary: {
      main: colors.secondary1[400],
      light: colors.secondary2[400],
      contrastText: colors.neutral[100],
      dark: colors.secondary1[600]
    },
    typography: {
      fontFamily: "" // TODO: choose between Open Sans, Roboto, Lato
    },
    error: {
      main: colors.error.main,
      dark: colors.error.dark
    },
    warning: {
      main: colors.warning.main,
      dark: colors.warning.dark
    },
    success: {
      main: colors.success.main,
      dark: colors.success.dark
    }
  }
})
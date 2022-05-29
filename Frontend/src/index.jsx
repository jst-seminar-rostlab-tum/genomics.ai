import React from 'react';
import ReactDOM from 'react-dom';
import './index.css';
import { CssBaseline } from '@mui/material';
import App from './App';
import { AuthProvider } from 'shared/context/authContext';
import { SubmissionProgressProvider } from 'shared/context/submissionProgressContext';

ReactDOM.render(
  <>
    <CssBaseline />
    <AuthProvider>
      <SubmissionProgressProvider>
        <App />
      </SubmissionProgressProvider>
    </AuthProvider>
  </>,
  document.getElementById('root'),
);

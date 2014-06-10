//
//  AppDelegate.m
//
//  Created by Ivan Andrus on 26/6/10.
//  Copyright 2010 __MyCompanyName__. All rights reserved.
//

#import "AppDelegate.h"
#import "AppController.h"
#import <WebKit/WebFrame.h>
#import <WebKit/WebView.h>
#import <Carbon/Carbon.h>
#import <AppKit/NSPasteboard.h>

@implementation AppDelegate

+ (void)initialize{
    // Make sure default are up to date
    NSUserDefaults *defaults = [NSUserDefaults standardUserDefaults];
    NSDictionary *factoryDefaults = [NSDictionary dictionaryWithContentsOfFile: [[NSBundle mainBundle] pathForResource:@"Defaults" ofType:@"plist"]];
    [defaults registerDefaults: factoryDefaults];
}

- (void)applicationWillFinishLaunching:(NSNotification *)aNotification{
    // This is early enough to show in the dock if we want to
    // http://www.cocoadev.com/index.pl?LSUIElement
    // http://codesorcery.net/2008/02/06/feature-requests-versus-the-right-way-to-do-it

    NSUserDefaults *defaults = [NSUserDefaults standardUserDefaults];
    BOOL isInDock = [defaults boolForKey:@"myShowInDock"];

    if ( isInDock ) {
        ProcessSerialNumber psn = { 0, kCurrentProcess };
        // display dock icon
        OSStatus returnCode = TransformProcessType(&psn, kProcessTransformToForegroundApplication);
        if( returnCode != 0 ) {
            // According to http://www.cocoadev.com/index.pl?TransformProcessType
            // TransformProcessType is available since 10.3, but doen't work for our case until 10.5
            NSLog(@"Could not show Sage.app in the dock. Error %d", (int)returnCode);
            // It's forbidden to showInDock since it doesn't work
            [defaults setBool:NO forKey:@"myShowInDock"];
            [defaults synchronize];

        } else {

            // enable menu bar
            SetSystemUIMode(kUIModeNormal, 0);
            // switch to Dock.app
            [[NSWorkspace sharedWorkspace] launchAppWithBundleIdentifier:@"com.apple.dock"
                                                                 options:NSWorkspaceLaunchDefault
                                          additionalEventParamDescriptor:nil
                                                        launchIdentifier:nil];
            // switch back
            [[NSApplication sharedApplication] activateIgnoringOtherApps:TRUE];
        }
    } else {
        // NSLog(@"Not showing in Dock");
    }

    // If we are using the system browser we don't need all of the menus
    // TODO: make this use menu titles not indexes
    if ( [defaults boolForKey:@"useSystemBrowser"] ) {
        [[NSApp mainMenu] removeItemAtIndex:6]; // Window menu
        [[NSApp mainMenu] removeItemAtIndex:2]; // Edit menu
    }

    [NSApp setServicesProvider:self];
}

- (void)applicationDidFinishLaunching:(NSNotification *)aNotification {
    // Register that we can open URLs
    NSAppleEventManager *em = [NSAppleEventManager sharedAppleEventManager];
    [em setEventHandler:self
            andSelector:@selector(getUrl:withReplyEvent:)
          forEventClass:kInternetEventClass
             andEventID:kAEGetURL];
}

- (BOOL)applicationShouldOpenUntitledFile:(NSApplication *)sender{
    return NO;
}

# pragma mark NSApplicationDelegate

- (NSString *)urlEncodeValue:(NSString *)str {
    NSString *result = (NSString *)
    CFURLCreateStringByAddingPercentEscapes(kCFAllocatorDefault,
                                            (CFStringRef)str,
                                            NULL,
                                            CFSTR(":/?#[]@!$&’()*+,;="),
                                            kCFStringEncodingUTF8);
    return [result autorelease];
}

// From here down are methods from NSApplicationDelegate, which probably do belong in another file.
// If/when this is done, I think you have to change the "File's Owner"'s delegate in IB
- (BOOL)application: (NSApplication * )theApplication openFile: (NSString * )filename{

    NSString *extension = [filename pathExtension];
    NSLog(@"Told to open %@ of type %@", filename, extension);

    // Handle the file based on extension
    if ( [extension isEqual:@"sage"] || [extension isEqual:@"py"] ) {
        // Run sage and python files in your terminal
        [appController sageTerminalRun:nil withArguments:[NSArray arrayWithObject:filename]];

    } else if ( [extension isEqual:@"sws"]
                || [extension isEqual:@"txt"]
                || [extension isEqual:@"zip"] )
    {

        // Browse to a url which will upload the file.
        // Perhaps we should have an option to delete the file when done...
        NSString* theURL = [NSString stringWithFormat:@"upload_worksheet?url=%@",
                           [[NSURL fileURLWithPath: filename] relativeString]];
        [appController browseLocalSageURL:theURL];

    } else if ( [extension isEqual:@"spkg"] ) {
        // Install the spkg
        [appController sageTerminalRun:@"i" withArguments:[NSArray arrayWithObject:filename]];

    } else if ( [extension isEqual:@"html"] || [extension isEqual:@"htm"] ) { // maybe others?

        NSError *outError = nil;
        id myDocument = [[NSDocumentController sharedDocumentController]
                         openUntitledDocumentAndDisplay:YES error:&outError];
        if ( myDocument == nil ) {
            [NSApp presentError:outError];
            NSLog(@"sageBrowser: Error creating document: %@", [outError localizedDescription]);
        } else {
            [[[myDocument webView] mainFrame]
             loadRequest:[NSURLRequest requestWithURL:
                          [NSURL URLWithString:
                           [NSString stringWithFormat:@"file://%@",
                            [filename stringByAddingPercentEscapesUsingEncoding:NSASCIIStringEncoding]]]]];
        }

    } else if ( [[NSFileManager defaultManager] isExecutableFileAtPath:[NSString stringWithFormat:@"%@/sage", filename]] ) {
        // Use this as the sage folder
        // NSFileManager *fileManager = [NSFileManager defaultManager];
        [[NSUserDefaults standardUserDefaults] setObject:[NSString stringWithFormat:@"%@/sage", filename]
                                                  forKey:@"SageBinary"];
        [appController setupPaths];

    } else {
        NSLog(@"I have no idea how I got a file of type %@.", extension);
        return NO;
    }
    return YES;
}

// http://stackoverflow.com/questions/49510/how-do-you-set-your-cocoa-application-as-the-default-web-browser
- (void)getUrl:(NSAppleEventDescriptor *)event withReplyEvent:(NSAppleEventDescriptor *)replyEvent{
    // Get the URL
    NSString *urlStr = [[event paramDescriptorForKeyword:keyDirectObject] stringValue];
    // Activate us
    [[NSApplication sharedApplication] activateIgnoringOtherApps:TRUE];

    if ( [[urlStr pathExtension] isEqual:@"spkg"] ) {
        // We can install spkg's from URLs
        [appController sageTerminalRun:@"i" withArguments:[NSArray arrayWithObject:urlStr]];

    } else if ( (  [[urlStr pathExtension] isEqual:@"sws"]
                || [[urlStr pathExtension] isEqual:@"txt"]
                || [[urlStr pathExtension] isEqual:@"zip"] )
               &&
               ! ( [[urlStr substringToIndex:16] isEqual:@"http://localhost"]
                || [[urlStr substringToIndex:17] isEqual:@"https://localhost"] ) )
    {

        // Browse to a url which will upload the file.
        // Perhaps we should have an option to delete the file when done...
        NSString* theURL = [NSString stringWithFormat:@"upload_worksheet?url=%@",
                            [[NSURL URLWithString: urlStr] relativeString]];
        [appController browseLocalSageURL:theURL];

    } else {

        // Open the url in a new window
        NSError *outError = nil;
        id myDocument = [[NSDocumentController sharedDocumentController]
                         openUntitledDocumentAndDisplay:YES error:&outError];
        if ( myDocument == nil ) {
            [NSApp presentError:outError];
            NSLog(@"sageBrowser: Error creating document: %@", [outError localizedDescription]);
        } else {
            [[[myDocument webView] mainFrame]
             loadRequest:[NSURLRequest requestWithURL:[NSURL URLWithString:urlStr]]];
        }

        // Check if the server has started
        // TODO: This detection will only work if we are SAGE_BROWSER (i.e. we are in the dock)
        NSArray *components = [urlStr componentsSeparatedByString:@"/?startup_token="];
        if ( [components count] > 1 ) {
            urlStr = [components objectAtIndex:0];
            components = [urlStr componentsSeparatedByString:@"localhost:"];
            if ( [components count] > 1 ) {
                const int port = (int)[[components objectAtIndex:1] floatValue];
                // We need to give it some time to load before we start loading queued things
                // which happens from serverStartedWithPort
                if ([[myDocument webView] respondsToSelector: @selector(isLoading)]) {
                    // block while the webview loads
                    while ([[myDocument webView] isLoading]) {
                        [[NSRunLoop currentRunLoop]
                         runMode:NSDefaultRunLoopMode
                         beforeDate:[NSDate distantFuture]];
                    }
                } else {
                    // Eyeball it...  This should only happen before 10.4.11
                    sleep(1);
                }
                [appController serverStartedWithPort:port];
            }
        }
    }
}


-(IBAction)openDocumentWithDialogBox:(id)sender{
    NSLog(@"openDocument:%@",sender);

    // Create the File Open Dialog class.
    NSOpenPanel* openDlg = [NSOpenPanel openPanel];
    [openDlg setCanChooseFiles:YES];
    [openDlg setCanChooseDirectories:NO];

    // Display the dialog.  If the OK button was pressed,
    // process the files.
    if ( [openDlg runModalForDirectory:nil file:nil] == NSOKButton )
    {
        // Get an array containing the full filenames of all
        // files and directories selected.
        NSArray* files = [openDlg filenames];
        for( int i = 0; i < [files count]; i++ )
        {
            NSString* fileName = [files objectAtIndex:i];
            [self application:nil openFile:fileName];
        }
    }
}


#pragma mark Services

- (NSArray*)extractFileURLs:(NSPasteboard *)pboard error:(NSString **)error {

    NSArray *classes = [NSArray arrayWithObject:[NSURL class]];

    // Only on 10.6
    if ( [pboard respondsToSelector:@selector(readObjectsForClasses:options:)] ) {

        NSArray* fileURLs = [pboard readObjectsForClasses:classes
                                                  options:nil];


        NSMutableArray* tmpArray = [NSMutableArray arrayWithCapacity:[fileURLs count]];
        NSEnumerator *enumerator = [fileURLs objectEnumerator];
        id key;
        while( (key = [enumerator nextObject]) ) {
            if ([key isFileURL]) {
                [tmpArray addObject:[key path]];
            }
        }
        NSArray *filePaths = [NSArray arrayWithArray:tmpArray];
        return filePaths;

    } else if ([[pboard types] containsObject:NSURLPboardType] ) {
        // On 10.5 I think you can only get one item no the paste board at a time, so this should be okay.
        // Ensure a URL on the pasteboard.
        NSURL *url = [NSURL URLFromPasteboard:pboard];
        if ([url isFileURL]) {
            return [NSArray arrayWithObject:[url path]];
        }
    }
    *error = NSLocalizedString(@"Error: couldn't execute file.",
                               @"pboard couldn't give URL.");
    return nil;
}


- (void)serviceTerminalRun:(NSPasteboard *)pboard userData:(NSString *)userData error:(NSString **)error {

    NSArray* files = [self extractFileURLs:pboard error:error];
    if (files) {
        [appController sageTerminalRun:userData withArguments:files];
    }
}

- (void)serviceExecute:(NSPasteboard *)pboard userData:(NSString *)userData error:(NSString **)error {
    NSLog(@"userData: %@",userData);
    NSArray *supportedTypes = [NSArray arrayWithObjects:NSStringPboardType, nil];
    NSString *bestType = [pboard availableTypeFromArray:supportedTypes];
    if ( bestType ) {
        NSString *str = [pboard stringForType:bestType];
        [appController sageTerminalRun: ([userData length] > 0 ? userData : nil)
                         withArguments:[NSArray
                                        arrayWithObjects:
                                        @"\n",
                                        str,
                                        nil]];
    }
}

- (void)serviceUnspkg:(NSPasteboard *)pboard userData:(NSString *)userData error:(NSString **)error {
    NSArray* files = [self extractFileURLs:pboard error:error];
    if (files) {

        NSEnumerator *enumerator = [files objectEnumerator];
        NSMutableArray *tasks = [NSMutableArray arrayWithCapacity:[files count]];

        // Launch the tasks
        NSTask *t;
        id file;
        while( (file = [enumerator nextObject]) ) {
            NSLog(@"unspking %@",file);
            [tasks addObject:[NSTask launchedTaskWithLaunchPath:@"/usr/bin/tar"
                                                      arguments:[NSArray arrayWithObjects:@"xjf", file,
                                                                 @"-C", [file stringByDeletingLastPathComponent],
                                                                 nil]]];
        }

        // Wait for all tasks to exit.
        enumerator = [tasks objectEnumerator];
        while( (t = [enumerator nextObject]) ) {
            [t waitUntilExit];
            NSLog(@"Exited with status %d", [t terminationStatus]);
        }
    }
}


- (NSString*)writeData:(NSData*)data toTempFileWithSuffix:(NSString*)suffix{

    NSString *tempFileTemplate = [NSTemporaryDirectory()
                                  stringByAppendingPathComponent:
                                  [NSString stringWithFormat:@"sage.XXXX%@", suffix]];
    const char *tempFileTemplateCString = [tempFileTemplate fileSystemRepresentation];
    char *tempFileNameCString = (char *)malloc(strlen(tempFileTemplateCString) + 1);
    strcpy(tempFileNameCString, tempFileTemplateCString);
    int fileDescriptor = mkstemps(tempFileNameCString, [suffix length]);
    if (fileDescriptor == -1) {
        free(tempFileNameCString);
        return nil;
    }

    NSString* tempFileName = [[NSFileManager defaultManager]
                              stringWithFileSystemRepresentation:tempFileNameCString
                              length:strlen(tempFileNameCString)];

    free(tempFileNameCString);

    NSFileHandle* tempFileHandle = [[NSFileHandle alloc]
                                    initWithFileDescriptor:fileDescriptor
                                    closeOnDealloc:YES];

    // write to file
    [tempFileHandle writeData:data];
    [tempFileHandle release];

    // Return path to the file
    return [NSString stringWithString:tempFileName];
}


- (void)serviceWorksheet:(NSPasteboard *)pboard userData:(NSString *)userData error:(NSString **)error {

    NSArray *supportedTypes = [NSArray arrayWithObjects:NSStringPboardType, nil];
    NSString *bestType = [pboard availableTypeFromArray:supportedTypes];
    if ( bestType ) {
        // Extract data from pasteboard
        NSString *str = [pboard stringForType:bestType];

        // Give it a title and surround with {{{ }}} if needed
        // If it contains {{{ we assume they know what they are doing and pass it straight through
        NSData* data;
        NSRange match = [str rangeOfString:@"{{{" options:NSLiteralSearch];
        if ( match.length == 0 ) {
            NSString* tmpStr = [NSString stringWithFormat:@"Untitled\n{{{\n%@\n}}}", str];
            data = [tmpStr dataUsingEncoding:NSUTF8StringEncoding];
        } else {
            data = [str dataUsingEncoding:NSUTF8StringEncoding];
        }

        // Write data to file
        NSString* tmpPath = [self writeData:data toTempFileWithSuffix:@".txt" ];

        // upload temp file
        if (tmpPath) {
            [self application:nil openFile:tmpPath];
        }
        // tmp file will get cleaned up in 3 days by the system...
    }
}


- (void)servicePreparse:(NSPasteboard *)pboard userData:(NSString *)userData error:(NSString **)error {
    NSArray *supportedTypes = [NSArray arrayWithObjects:NSStringPboardType, nil];
    NSString *bestType = [pboard availableTypeFromArray:supportedTypes];
    if ( bestType ) {
        // Extract data from pasteboard
        NSString *str = [pboard stringForType:bestType];

        // Write data to a temporary file
        NSData* data = [str dataUsingEncoding:NSUTF8StringEncoding];
        NSString* tmpPath = [self writeData:data toTempFileWithSuffix:@".sage" ];

        // preparse temp file
        if (tmpPath) {
            NSTask* preparser = [NSTask launchedTaskWithLaunchPath:[appController sageBinary]
                                                         arguments:[NSArray arrayWithObjects:@"--preparse",
                                                                    tmpPath,
                                                                    nil]];
            [preparser waitUntilExit];
            if ( [preparser terminationStatus] == 0 ) {
                NSString* preparsedFile = [[tmpPath stringByDeletingPathExtension]
                                           stringByAppendingPathExtension:@"py"];
                NSError* err;
                NSString* contents = [NSString stringWithContentsOfFile:preparsedFile
                                                               encoding:NSUTF8StringEncoding
                                                                  error:&err];
                if (contents) {
                    // The first line contains an autogenerated message that we don't want
                    NSRange firstLine = [contents rangeOfString:@"\n"];
                    NSRange autoGenerated = [contents rangeOfString:@"This file was *autogenerated* from the file"];
                    [pboard declareTypes:supportedTypes owner:nil];

                    if ( autoGenerated.length > 0 && firstLine.length > 0 && autoGenerated.location < firstLine.location ) {
                        [pboard setString:[contents substringFromIndex:firstLine.location+firstLine.length] forType:NSStringPboardType];
                    } else {
                        [pboard setString:contents forType:NSStringPboardType];
                    }
                } else if (error) {
                    *error = [err localizedDescription];
                }
            } else if (error) {
                *error = @"Could not preparse file";
            }
        } else if (error) {
            *error = @"Could not create temp file";
        }
        // tmp file will get cleaned up in 3 days by the system...
    }
}


- (void)serviceEval:(NSPasteboard *)pboard userData:(NSString *)userData error:(NSString **)error {

    NSArray *supportedTypes = [NSArray arrayWithObjects:NSStringPboardType, nil];
    NSString *bestType = [pboard availableTypeFromArray:supportedTypes];
    if ( bestType ) {
        NSString *str = [pboard stringForType:bestType];

        NSPipe* sageOut = [NSPipe pipe];
        NSTask* task = [[NSTask alloc] init];
        [task setLaunchPath:[appController sageBinary]];
        [task setArguments:[NSArray arrayWithObjects:@"-c", str, nil]];
        [task setStandardOutput:sageOut];

        [task launch];

        // It it writes more than 8k we may be in trouble...
        [task waitUntilExit];

        if ([task terminationStatus] == 0 ) {

            NSData *data = [[sageOut fileHandleForReading] readDataToEndOfFile];
            if ( [data length] == 0 ) {
                data = [[NSString
                         stringWithFormat:@"%@\nThe code produced no output.  Please use print or similar to see the result.\n",
                         str] dataUsingEncoding:NSUTF8StringEncoding];
            }
            [pboard declareTypes:supportedTypes owner:nil];
            [pboard setData:data forType:NSStringPboardType];

        } else if (error) {
            *error = @"Could not evaluate code.";
        }

        [task release];
    }
}


@end
